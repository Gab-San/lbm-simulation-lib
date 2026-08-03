#ifndef __LBM_SIM_SOLVER_SOLVER_2D
#define __LBM_SIM_SOLVER_SOLVER_2D

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/backend/metadata.hpp"
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

// COLLISION DETECTION LIB
#include "collision-detection/core/operators.hpp"
#include "collision-detection/core/types.hpp"

// C++ STANDARD LIB
#include <cmath>
#include <cstdint>
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
  MPISolver2D(const unsigned int num_iters_, const unsigned int num_frames_,
              const Solid::types::boundary_mask_t &boundary_mask_)
      : Base(num_iters_, num_frames_, boundary_mask_) {};

  ~MPISolver2D() = default;

  // TODO: adapt to initialization
  void init_equilibrium(const Grid<2> &grid,
                        std::vector<double> &part_stream) const override {
#pragma omp parallel for shared(grid, part_stream) schedule(static) collapse(2)
    for (unsigned int y = 0; y < grid.size.y; ++y) {
      for (unsigned int x = 0; x < grid.size.x; ++x) {
        const CollisionDetection::types::Coordinate<2> p(x, y);

        double r = grid.rho[grid.scalar_index(p)];
        const CollisionDetection::utils::Vector<double, 2> u =
            grid.u[grid.scalar_index(p)];
        const double u_sq = CollisionDetection::utils::dot(u, u);

#pragma omp simd
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          double cidotu = CollisionDetection::utils::dot(D2Q9::dir[i], u);

          part_stream[grid.field_index(p, i, D2Q9::ndir)] =
              D2Q9::wi[i] * r *
              (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
        }
      }
    }
  };

  void solve(Grid<2> &grid, const Params<2, cm_t> &params_,
             std::vector<double> &ffrom,
             std::vector<double> &fto) const override {
    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);

    for (unsigned int iter = 0; iter < this->niters; iter++) {
      bool save = (iter % this->nskips == 0);
      update_stream_collide(grid, cs, ffrom, fto, save);
      std::swap(ffrom, fto);
      if (save)
        write_norms(grid);
    }
  }

private:
  // FIXME: Execution context ??
  void update_stream_collide(
      Grid<2> &grid, const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> &cs,
      const std::vector<double> &ffrom, std::vector<double> &fto, bool save,
      const ExecutionContext<OPEN_MP> &context =
          ExecutionContext<OPEN_MP>{}) const {

    // STREAMING
#pragma omp parallel for shared(ffrom, fto, cs, grid, save) schedule(static)   \
    collapse(2)
    for (std::size_t y = 0; y < grid.size.y; ++y) {
      for (std::size_t x = 0; x < grid.size.x; ++x) {
        std::array<double, D2Q9::ndir> fp;
        const CollisionDetection::utils::Point<int, 2> p(x, y);

#pragma omp unroll full
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          const CollisionDetection::types::Coordinate<2> src = p - D2Q9::dir[i];

          // stream only if the source node is internal
          if (!grid.contains(src)) {
            fp[i] = ffrom[grid.field_index(p, i, D2Q9::ndir)];
          } else {
            fp[i] = ffrom[grid.field_index(src, i, D2Q9::ndir)];
          }
        }

        // FULLWAY-BOUNCE-BACK

        double r = 0.0;
#pragma omp unroll full
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          // calculate local rho before collision
          r += fp[i];
        }

        apply_boundary_conditions(p, r, cs.params.init_vel, fp, grid, context);

        // COMPUTE MACROSCOPIC VARIABLES

        // rho = sum_i fi
        // rho*u = sum_i fi * ci

        r = 0.0;
        CollisionDetection::utils::Vector<double, 2> u(0, 0);

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
          const unsigned int s_idx = grid.scalar_index(p);
          grid.rho[s_idx] = r;
          grid.u[s_idx] = u;
        }

        // APPLY COLLISION
        cs.apply(p, u, r, fp, grid);

        // COPY LOCAL DENSITY TO GRID
        for (auto i = 0; i < D2Q9::ndir; i++) {
          fto[grid.field_index(p, i, D2Q9::ndir)] = fp[i];
        }
      }
    }
  };

  void apply_boundary_conditions(
      const CollisionDetection::types::Coordinate<2> p, const double &localrho,
      const CollisionDetection::utils::Vector<double, 2> u0,
      std::array<double, D2Q9::ndir> &fp, const Grid<2> &grid,
      const ExecutionContext<OPEN_MP> &context =
          ExecutionContext<OPEN_MP>{}) const {
    (void)context;

    Solid::types::boundary_t b = this->strt[grid.scalar_index(p)];

#pragma omp unroll full
    for (std::size_t diridx = 0; diridx < D2Q9::ndir; diridx++) {
      if (!grid.contains(p + D2Q9::dir[diridx])) {
        switch (b) {
        case Solid::BB_RIGID_WALL:
          Solid::apply_bb_rigid_wall<D2Q9>(fp, diridx);
          break;
        case Solid::BB_MOVING_WALL:
          Solid::apply_bb_moving_wall<2, D2Q9>(localrho, u0, fp, diridx);
          break;
        default:
          break;
        }
      }
    }
  }

  // Il calcolo delle norme resta interno al solver (come da checklist);
  // cambia solo la destinazione dei dati calcolati: invece di scrivere
  // direttamente su un AsyncBinaryWriter posseduto dal solver, notifica
  // i listener registrati tramite DataObservable (pattern Observer/Listener).
  void write_norms(const Grid<2> &grid) const override {
    using CollisionDetection::utils::dot;
    std::vector<float> vsq(grid.getArea());

#pragma omp parallel for collapse(2)
    for (unsigned int y = 0; y < grid.size.y; ++y) {
      for (unsigned int x = 0; x < grid.size.x; ++x) {
        const CollisionDetection::types::Coordinate<2> p(x, y);
        const auto &vel = grid.u[grid.scalar_index(p)];

        vsq[grid.scalar_index(p)] =
            static_cast<float>(std::sqrt(dot(vel, vel)));
      }
    }

    std::vector<char> buf(vsq.size() * sizeof(float));
    std::memcpy(buf.data(), vsq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }
}; // class MPISolver2D

} // namespace lbm
#endif // __LBM_SIM_SOLVER_SOLVER_2D
