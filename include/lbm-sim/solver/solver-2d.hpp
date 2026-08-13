#ifndef __LBM_SIM_SOLVER_SOLVER_2D
#define __LBM_SIM_SOLVER_SOLVER_2D

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/backend/metadata.hpp"
#include "lbm-sim/boundaries.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"

#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

// C++ STANDARD LIB
#include <cmath>
#include <cstring>
#include <vector>

#include <assert.h>

// OMP LIB
#include <omp.h>

namespace lbm {
template <enum CollisionModel cm_t>
class MPISolver2D : public SolverBase2D<cm_t, ExecutionBackend::OPEN_MP> {
  using Base = SolverBase2D<cm_t, ExecutionBackend::OPEN_MP>;

public:
  MPISolver2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {};

  ~MPISolver2D() = default;

  // TODO: adapt to initialization
  void init_equilibrium(const Lattice<2> &lattice,
                        std::vector<double> &part_stream) const override {
    using utils::ops::dot;
#pragma omp parallel for shared(lattice, part_stream) schedule(static)         \
    collapse(2)
    for (unsigned int y = 0; y < lattice.grid.size.y; ++y) {
      for (unsigned int x = 0; x < lattice.grid.size.x; ++x) {
        const types::Coordinate<2> p(x, y);

        double r = lattice.rho[lattice.grid.scalar_index(p)];
        const utils::Vector<double, 2> u =
            lattice.u[lattice.grid.scalar_index(p)];
        const double u_sq = dot(u, u);

#pragma omp simd
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          double cidotu = dot(D2Q9::dir[i], u);

          part_stream[lattice.grid.field_index(p, i, D2Q9::ndir)] =
              D2Q9::wi[i] * r *
              (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
        }
      }
    }
  };

  void solve(Lattice<2> &lattice, const Params<2, cm_t> &params_,
             std::vector<double> &ffrom,
             std::vector<double> &fto) const override {
    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);
    // partono iterazioni
    for (unsigned int iter = 0; iter < this->niters; iter++) {
      bool save = (iter % this->nskips == 0);
      update_stream_collide(lattice, cs, ffrom, fto, save);
      std::swap(ffrom, fto);
      if (save)
        write_norms(lattice);
    }
  }

private:
  // FIXME: Execution context ??
  void update_stream_collide(
      Lattice<2> &lattice, const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> &cs,
      const std::vector<double> &ffrom, std::vector<double> &fto, bool save,
      const ExecutionContext<OPEN_MP> &context =
          ExecutionContext<OPEN_MP>{}) const {

    // STREAMING + HALFWAY COLLISION
#pragma omp parallel for shared(ffrom, fto, cs, lattice, save)                 \
    schedule(static) collapse(2)
    for (std::size_t y = 0; y < lattice.grid.size.y; ++y) {
      for (std::size_t x = 0; x < lattice.grid.size.x; ++x) {
        std::array<double, D2Q9::ndir> fp;
        const utils::Point<int, 2> p(x, y);

        double r_wall = 0.0;
#pragma omp unroll full
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          // calculate local rho on wall before boundary conditions
          r_wall += ffrom[lattice.grid.field_index(p, i, D2Q9::ndir)];
        }

#pragma omp unroll full
        for (unsigned int diridx = 0; diridx < D2Q9::ndir; ++diridx) {
          const types::Coordinate<2> src = p - D2Q9::dir[diridx];

          if (!lattice.grid.contains(src)) {
            // if source node is external it is on a boundary node
            apply_boundary_conditions(lattice.boundary_mask, fp, ffrom, diridx,
                                      lattice.grid, p, r_wall,
                                      cs.params.init_vel, context);
          } else {
            // if source node is internal stream it.
            fp[diridx] =
                ffrom[lattice.grid.field_index(src, diridx, D2Q9::ndir)];
          }
        }

        // COMPUTE MACROSCOPIC VARIABLES

        // rho = sum_i fi
        // rho*u = sum_i fi * ci

        double r = 0.0;
        utils::Vector<double, 2> u(0, 0);

#pragma omp unroll full
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {

          // calculate macroscopic variables
          r += fp[i];
          u += D2Q9::dir[i] * fp[i];
        }

        // u = (sum_i fi * ci) / rho
        u /= r;

        // STORE computed macroscopic values
        if (save) {
          const unsigned int s_idx = lattice.grid.scalar_index(p);
          lattice.rho[s_idx] = r;
          lattice.u[s_idx] = u;
        }

        // APPLY COLLISION
        cs.apply(p, u, r, fp, lattice.grid);

        // COPY LOCAL DENSITY TO GRID
#pragma omp simd
        for (auto i = 0; i < D2Q9::ndir; i++) {
          fto[lattice.grid.field_index(p, i, D2Q9::ndir)] = fp[i];
        }
      }
    }
  };

  void apply_boundary_conditions(const types::boundary_mask_t &boundary_mask,
                                 std::array<double, D2Q9::ndir> &fp,
                                 const std::vector<double> &ffrom,
                                 const std::size_t diridx, const Grid<2> &grid,
                                 const types::Coordinate<2> p,
                                 const double &localrho,
                                 const utils::Vector<double, 2> u0,
                                 const ExecutionContext<OPEN_MP> &context =
                                     ExecutionContext<OPEN_MP>{}) const {
    (void)context;

    types::boundary_t b = boundary_mask[grid.scalar_index(p)];

    switch (b) {
    case Solid::BB_RIGID_WALL:
      Solid::apply_bb_rigid_wall<2, D2Q9>(fp, ffrom, diridx, grid, p);
      break;
    case Solid::BB_MOVING_WALL:
      Solid::apply_bb_moving_wall<2, D2Q9>(fp, ffrom, diridx, grid, p, localrho,
                                           u0);
      break;
    case Solid::PERIODIC:
      Solid::apply_continue<2, D2Q9>(fp, ffrom, diridx, grid, p, localrho);
      break;
    default:
      break;
    }
  }

  // Il calcolo delle norme resta interno al solver (come da checklist);
  // cambia solo la destinazione dei dati calcolati: invece di scrivere
  // direttamente su un AsyncBinaryWriter posseduto dal solver, notifica
  // i listener registrati tramite DataObservable (pattern Observer/Listener).
  void write_norms(const Lattice<2> &lattice) const override {
    using utils::ops::dot;
    std::vector<float> vsq(lattice.grid.getArea());

#pragma omp parallel for collapse(2)
    for (unsigned int y = 0; y < lattice.grid.size.y; ++y) {
      for (unsigned int x = 0; x < lattice.grid.size.x; ++x) {
        const types::Coordinate<2> p(x, y);
        const auto &vel = lattice.u[lattice.grid.scalar_index(p)];

        vsq[lattice.grid.scalar_index(p)] =
            static_cast<float>(std::sqrt(dot(vel, vel)));
      }
    }

    std::vector<char> buf(vsq.size() * sizeof(float));
    std::memcpy(buf.data(), vsq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }

}; // class MPISolver2D

} // namespace lbm
/*
namespace poiseuille {

// Calcola u_max dato il gradiente di pressione imposto e la viscosità
inline double compute_u_max(double rho_in, double rho_out, double Lx,
                             double H, double tau, double cs2 = 1.0/3.0) {
    double dp_dx = (rho_in - rho_out) * cs2 / Lx;
    double nu = cs2 * (tau - 0.5);
    return dp_dx * H * H / (8.0 * nu);
}

// Profilo analitico di velocità in funzione di y
inline double analytic_ux(double y, double H, double u_max) {
    double y_c = y - H * 0.5;
    return u_max * (1.0 - (y_c * y_c) / (H * 0.5 * H * 0.5));
}

// Inizializzazione IC coerente col profilo (riduce drasticamente
// il transiente iniziale, come discusso)
template <typename GridT>
void initialize(GridT& grid, double rho0, double u_max, double H) {
    for (std::size_t y = 0; y < grid.size.y; ++y) {
        double ux = analytic_ux(static_cast<double>(y), H, u_max);
        for (std::size_t x = 0; x < grid.size.x; ++x) {
            auto p = types::Coordinate<2>(x, y);
            auto feq = bc::equilibrium(rho0, ux, 0.0);
            // scrivi feq nelle f del nodo (x,y)
            for (int i = 0; i < 9; ++i) {
                grid.f[grid.field_index(p, i, 9)] = feq[i];
            }
            grid.rho[grid.scalar_index(p)] = rho0;
            grid.u[grid.scalar_index(p)] = {ux, 0.0};
        }
    }
}

} // namespace poiseuille
 */
#endif // __LBM_SIM_SOLVER_SOLVER_2D
