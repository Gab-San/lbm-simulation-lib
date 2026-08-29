#ifndef __LBM_SIM_SOLVER_SOLVER_2D
#define __LBM_SIM_SOLVER_SOLVER_2D

#include "lbm-sim/boundaries/boundary-conditions.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/omp/annotations.hpp"
#include "lbm-sim/omp/iteration.hpp"
#include "lbm-sim/solver/solver-base.hpp"
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

namespace lbm {

template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
class OpenMPSolver
    : public SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::OPEN_MP> {
  using Base = SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::OPEN_MP>;

public:
  OpenMPSolver(const unsigned int iters_, const unsigned int frames_)
      : Base(iters_, frames_) {};

  ~OpenMPSolver() = default;

  void solve(Lattice<dim> &lattice,
             const CollisionParams<dim, cm_t> &params_) const override {
    quill::Logger *solver_logger = logging::create_or_get_logger("solver");

    std::vector<double> ffrom(lattice.grid.getArea() * VelocitySet::ndir, 0.0);
    std::vector<double> fto(lattice.grid.getArea() * VelocitySet::ndir, 0.0);
    std::vector<float> usq(lattice.grid.getArea());

    const CollisionStrategy<dim, VelocitySet, cm_t> cs(params_);

    init_equilibrium(ffrom, lattice);

    LOG_DEBUG(solver_logger, "Equilibrium Initialized...");

    LOG_INFO(solver_logger, "System has {} logical processors.",
             omp_get_num_procs());
    LOG_INFO(solver_logger, "The parallel section will run on {} threads.",
             omp_get_max_threads());

    for (std::size_t iter = 0; iter < this->niters; iter++) {
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
  inline void init_equilibrium(std::vector<double> &part_stream,
                               const Lattice<dim> &lattice) const {
    const auto ext = lattice.grid.extents();
    const auto area = static_cast<std::ptrdiff_t>(lattice.grid.getArea());

    using utils::ops::dot;

#pragma omp parallel for shared(lattice, part_stream) schedule(runtime)
    for (int cell = 0; cell < area; cell++) {

      const types::Coordinate<dim> p = iteration::unflatten<dim>(cell, ext);

      double r = lattice.rho[cell];
      const utils::Vector<double, dim> u = lattice.u[cell];
      const double u_sq = dot(u, u);

#pragma omp simd
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); ++diridx) {
        double cidotu = dot(VelocitySet::dir[diridx], u);

        part_stream[lattice.grid.field_index(p, diridx, VelocitySet::ndir)] =
            VelocitySet::wi[diridx] * r *
            (1.0 + numbers::invcs_2 * cidotu + 4.5 * cidotu * cidotu -
             1.5 * u_sq);
      }
    }
  }

  void
  update_stream_collide(const std::vector<double> &ffrom,
                        std::vector<double> &fto, std::vector<float> &usq,
                        Lattice<dim> &lattice,
                        const CollisionStrategy<dim, VelocitySet, cm_t> &cs,
                        const bool store_macroscopic) const {
    const auto ext = lattice.grid.extents();
    const auto area = static_cast<std::ptrdiff_t>(lattice.grid.getArea());

#pragma omp parallel for shared(ffrom, fto, cs, lattice, store_macroscopic)    \
    schedule(static)
    for (int cell = 0; cell < area; ++cell) {
      std::array<double, VelocitySet::ndir> fp;
      const types::Coordinate<dim> p = iteration::unflatten<dim>(cell, ext);

      // Skip solid nodes.
      //
      // Test solid_mask, NEVER a BC type: a fluid node sitting on a domain
      // edge carries a face BC and must not be skipped.
      //
      // No solid node's populations are ever read under this scheme --
      // bounce-back reads the fluid node p, the pressure rescale reads a node
      // already confirmed fluid, and plain streaming only runs when src is
      // fluid. So no surface/interior distinction is needed here.
      //
      // CAVEAT: interpolated bounce-back (Bouzidi, Filippova-Hanel) and
      // Guo-style extrapolation DO read the solid node. Adding either brings
      // the surface solid-node distinction back and this line has to change
      // with it.
      if (lattice.solid_mask[cell] != types::FLUID)
        continue;

      double r_wall = 0.0;
      UNROLL_FULL
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); ++diridx) {
        // calculate local rho on wall before boundary conddiridxtions
        r_wall += ffrom[lattice.grid.field_index(p, diridx, VelocitySet::ndir)];
      }

      // STREAMING + HALFWAY COLLISION
      UNROLL_FULL
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); ++diridx) {
        // One resolve_link per direction: domain faces, periodic wrap and
        // immersed obstacles are all decided in there, per link, not per node.
        const auto link = Solid::resolve_link<dim>(
            lattice.grid, lattice.domain_bc, lattice.solid_mask.data(),
            lattice.obstacles.data(), p, VelocitySet::dir[diridx]);

        if (link.bc == Solid::NONE) {
          // source node is fluid and in range: plain streaming.
          fp[diridx] = ffrom[lattice.grid.field_index(link.src, diridx,
                                                      VelocitySet::ndir)];
        } else {
          Solid::apply_boundary_condition<dim, VelocitySet>(
              fp.data(), ffrom.data(), diridx, lattice.grid, link,
              lattice.obstacles.data(), lattice.rho.data(), lattice.u.data(), p,
              r_wall, cs.params.init_vel, lattice.pin, lattice.pout);
        }
      }

      // COMPUTE MACROSCOPIC VARIABLES

      // rho = sum_i fi
      // rho*u = sum_i fi * ci

      double r = 0.0;
      utils::Vector<double, dim> u;

      UNROLL_FULL
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); ++diridx) {
        r += fp[diridx];
        u += VelocitySet::dir[diridx] * fp[diridx];
      }

      // u = (sum_i fi * ci) / rho
      u /= r;

      // STORE computed macroscopic values
      if (store_macroscopic ||
          Solid::on_pressure_face(lattice.grid, lattice.domain_bc, p)) {
        lattice.rho[cell] = r;
        lattice.u[cell] = u;
      }

      if (store_macroscopic) {
        usq[cell] = static_cast<float>(std::sqrt(utils::ops::dot(u, u)));
      }

      // APPLY COLLISION
      cs.apply(fp.data(), p, r, u);

      // COPY LOCAL DENSITY TO GRID
#pragma omp simd
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); diridx++) {
        fto[lattice.grid.field_index(p, diridx, VelocitySet::ndir)] =
            fp[diridx];
      }
    }
  };

  inline void write_norms(const std::vector<float> &usq) const {
    std::vector<char> buf(usq.size() * sizeof(float));
    std::memcpy(buf.data(), usq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm
#endif // __LBM_SIM_SOLVER_SOLVER_2D
