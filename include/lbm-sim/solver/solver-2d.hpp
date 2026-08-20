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

#include "lbm/logging.hpp"

#include "quill/LogMacros.h"

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

// FIXME: CHANGE NAME TO OpenMPSolver2D
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

    quill::Logger *solver_logger = logging::create_or_get_logger("solver");

    std::vector<double> ffrom(lattice.grid.getArea() * D2Q9::ndir, 0.0);
    std::vector<double> fto(lattice.grid.getArea() * D2Q9::ndir, 0.0);
    std::vector<float> usq(lattice.grid.getArea());

    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);

    init_equilibrium(ffrom, lattice);

    LOG_DEBUG(solver_logger, "Equilibrium Initialized...");

    LOG_INFO(solver_logger, "System has {} logical processors.",
             omp_get_num_procs());
    LOG_INFO(solver_logger, "System can work with up to {} threads.",
             omp_get_max_threads());

    for (unsigned int iter = 0; iter < this->niters; iter++) {
      bool save = iter % this->nskips == 0;
      update_stream_collide(ffrom, fto, usq, lattice, cs, save);
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
      const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> &cs, bool save) const {

#pragma omp parallel for shared(ffrom, fto, cs, lattice, save)                 \
    schedule(static) collapse(2)
    //for (std::size_t y = 0; y < lattice.grid.size.y; ++y) {
    //  for (std::size_t x = 0; x < lattice.grid.size.x; ++x) {
    // meglio usare auto y = 0; y < lattice.grid.size.y; ++y; ++y
    for (auto y = 0; y < lattice.grid.size.y; ++y) {
      for (auto x = 0; x < lattice.grid.size.x; ++x) {
        std::array<double, D2Q9::ndir> fp;
        const utils::Point<int, 2> p(x, y);

        LOG_TRACE_L3(logging::create_or_get_logger("openmp-iteration"),
                     "Running simulation on {}", p);

        double r_wall = 0.0;
// #pragma omp unroll full
   UNROLL_FULL
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          // calculate local rho on wall before boundary conditions
          r_wall += ffrom[lattice.grid.field_index(p, i, D2Q9::ndir)];
        }

        // STREAMING + HALFWAY COLLISION
// #pragma omp unroll full
    UNROLL_FULL 
        for (auto diridx = 0; diridx < D2Q9::ndir; ++diridx) {
          const types::Coordinate<2> src = p - D2Q9::dir[diridx];

          if (!lattice.grid.contains(src)) {
            // if source node is external it is on a boundary node
            apply_boundary_conditions(lattice.boundary_mask, fp, ffrom, diridx,
                                      lattice, p, r_wall, cs.params.init_vel);
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

// #pragma omp unroll full
    UNROLL_FULL
        for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
          // calculate macroscopic variables
          r += fp[i];
          u += D2Q9::dir[i] * fp[i];
        }

        // u = (sum_i fi * ci) / rho
        u /= r;

        // STORE computed macroscopic values
        const auto s_idx = lattice.grid.scalar_index(p);
        if (save ||
            lattice.boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_INLET ||
            lattice.boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_OUTLET) {
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

  inline void apply_boundary_conditions(
      const types::boundary_mask_t &boundary_mask,
      std::array<double, D2Q9::ndir> &fp, const std::vector<double> &ffrom,
      const std::size_t diridx, const Lattice<2> &lattice,
      const types::Coordinate<2> p, const double &localrho,
      const utils::Vector<double, 2> u0,
      const ExecutionContext<OPEN_MP> &context =
          ExecutionContext<OPEN_MP>{}) const {
    (void)context;

    types::boundary_t b = boundary_mask[lattice.grid.scalar_index(p)];

    switch (b) {
    case Solid::BB_RIGID_WALL:
      Solid::apply_bb_rigid_wall<2, D2Q9>(fp, ffrom, diridx, lattice.grid, p);
      break;
    case Solid::BB_MOVING_WALL:
      Solid::apply_bb_moving_wall<2, D2Q9>(fp, ffrom, diridx, lattice.grid, p,
                                           localrho, u0);
      break;
    case Solid::PERIODIC:
      Solid::apply_periodic<2, D2Q9>(fp, ffrom, diridx, lattice.grid, p);
      break;
    case Solid::PRESSURE_PERIODIC_INLET:
      Solid::apply_periodic_with_pressure_variation<2, D2Q9>(
          fp.data(), ffrom.data(), diridx, lattice.grid, p, lattice.rho.data(),
          lattice.pin, lattice.u.data());
      break;
    case Solid::PRESSURE_PERIODIC_OUTLET:
      Solid::apply_periodic_with_pressure_variation<2, D2Q9>(
          fp.data(), ffrom.data(), diridx, lattice.grid, p, lattice.rho.data(),
          lattice.pout, lattice.u.data());
      break;
    default:
      // FIXME: Should this throw an error?
      break;
    }
  }

  inline void write_norms(const std::vector<float> &usq) const {
    std::vector<char> buf(usq.size() * sizeof(float));
    std::memcpy(buf.data(), usq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }

}; // class MPISolver2D

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_2D
