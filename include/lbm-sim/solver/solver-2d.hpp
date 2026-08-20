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
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

// OMP LIB
#include <omp.h>

// Manage loop unrolling pragmas for different compilers
#pragma once
#if defined(__clang__)
    #define UNROLL_FULL _Pragma("clang loop unroll(full)")
#elif defined(__GNUC__)
    #define UNROLL_FULL _Pragma("GCC unroll 65534")
#else
    #define UNROLL_FULL
#endif

namespace lbm {

template <enum CollisionModel cm_t>
class MPISolver2D : public SolverBase2D<cm_t, ExecutionBackend::OPEN_MP> {
  using Base = SolverBase2D<cm_t, ExecutionBackend::OPEN_MP>;
  static constexpr unsigned short int DIM = 2;

public:
  MPISolver2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {};

  virtual ~MPISolver2D() = default;

  void solve(Lattice<2> &lattice,
             const Params<2, cm_t> &params_) const override {
    std::vector<double> ffrom(lattice.grid.getArea() * D2Q9::ndir, 0.0);

    init_equilibrium(ffrom, lattice);
    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);

    std::vector<double> fto(lattice.grid.getArea() * D2Q9::ndir, 0.0);
    std::vector<float> usq(lattice.grid.getArea());

    for (unsigned int iter = 0; iter < this->niters; iter++) {
      const bool save = iter % this->nskips == 0;
      const bool store_macroscopic = save || (iter + 1 == this->niters);
      update_stream_collide(ffrom, fto, usq, lattice, cs, store_macroscopic);
      std::swap(ffrom, fto);
      if (save) {
        write_norms(usq);
      }
    }
  }

private:
  // TODO: adapt to initialization
  inline void init_equilibrium(std::vector<double> &part_stream,
                               const Lattice<2> &lattice) const {
    using utils::ops::dot;
#pragma omp parallel for shared(lattice, part_stream) schedule(static)         \
    collapse(2)
    for (auto y = 0; y < lattice.grid.size.y; ++y) {
      for (auto x = 0; x < lattice.grid.size.x; ++x) {
        const types::Coordinate<2> p(x, y);

        double r = lattice.rho[lattice.grid.scalar_index(p)];
        const utils::Vector<double, 2> u =
            lattice.u[lattice.grid.scalar_index(p)];
        const double u_sq = dot(u, u);

#pragma omp simd
        for (auto i = 0; i < D2Q9::ndir; ++i) {
          double cidotu = dot(D2Q9::dir[i], u);

          part_stream[lattice.grid.field_index(p, i, D2Q9::ndir)] =
              D2Q9::wi[i] * r *
              (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
        }
      }
    }
  };

  // FIXME: Execution context ??
  void update_stream_collide(
      const std::vector<double> &ffrom, std::vector<double> &fto,
      std::vector<float> &usq, Lattice<2> &lattice,
      const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> &cs,
      const bool store_macroscopic) const {

#pragma omp parallel for shared(ffrom, fto, cs, lattice, store_macroscopic)    \
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

        // STREAMING + HALFWAY COLLISION
#pragma omp unroll full
        for (auto diridx = 0; diridx < D2Q9::ndir; ++diridx) {
          const types::Coordinate<2> src = p - D2Q9::dir[diridx];

          if (!lattice.grid.contains(src)) {
            // if source node is external it is on a boundary node
            Solid::apply_boundary_condition<2, D2Q9>(
                lattice.boundary_mask.data(), fp.data(), ffrom.data(), diridx,
                lattice.grid, p, r_wall, cs.params.init_vel);
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
          r += fp[i];
          u += D2Q9::dir[i] * fp[i];
        }

        // u = (sum_i fi * ci) / rho
        u /= r;

        // STORE computed macroscopic values
        if (store_macroscopic) {
          const unsigned int s_idx = lattice.grid.scalar_index(p);
          lattice.rho[s_idx] = r;
          lattice.u[s_idx] = u;
          usq[s_idx] = static_cast<float>(std::sqrt(utils::ops::dot(u, u)));
        }

        // APPLY COLLISION
        cs.apply(fp, p, r, u);

        // COPY LOCAL DENSITY TO GRID
#pragma omp simd
        for (auto i = 0; i < D2Q9::ndir; i++) {
          fto[lattice.grid.field_index(p, i, D2Q9::ndir)] = fp[i];
        }
      }
    }
  };


  inline void write_norms(const std::vector<float> &usq) const {
    std::vector<char> buf(usq.size() * sizeof(float));
    std::memcpy(buf.data(), usq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }

}; // class MPISolver2D

} // namespace lbm

/*
// Calcola u_max dato il gradiente di pressione imposto e la viscosità
inline double compute_u_max(double rho_in, double rho_out, double Lx, double H,
                            double tau, double cs2 = 1.0 / 3.0) {
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
void initialize(GridT &grid, double rho0, double u_max, double H) {
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
