/**
 * @file openmp-solver.hpp
 * @brief OpenMPSolver: the multi-threaded CPU backend.
 *
 * The whole time step is one fused pass, update_stream_collide(), over two
 * population buffers that are swapped at the end of the iteration:
 * populations are read from one and written to the other, so there is no
 * aliasing, no halo copy and no second sweep.
 *
 * Parallelism is a single `omp parallel for` over the flattened cell index,
 * `schedule(static)`. Every node is independent within a step -- a node only
 * ever reads its neighbours' *previous* populations -- so no synchronisation
 * is needed inside the pass.
 */

#ifndef __LBM_SIM_SOLVER_SOLVER_2D
#define __LBM_SIM_SOLVER_SOLVER_2D

#include "lbm-sim/backend/omp/annotations.hpp"
#include "lbm-sim/backend/omp/iteration.hpp"
#include "lbm-sim/backend/properties.hpp"
#include "lbm-sim/boundaries/boundary-conditions.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/formatting.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/profiling.hpp"
#include "lbm-sim/solver/solver-base.hpp"
#include "lbm/format/csv-writer.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

// OMP LIB
#include <omp.h>

namespace lbm {

/**
 * @brief Multi-threaded CPU solver.
 *
 * @tparam dim         Spatial dimension (2 or 3).
 * @tparam VelocitySet Discrete velocity set (D2Q9, D3Q19, D3Q27).
 * @tparam cm_t        Collision model.
 *
 * @see SolverBase for the iteration and frame bookkeeping, and
 *      CUDASolver for the GPU counterpart.
 */
template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
class OpenMPSolver
    : public SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::OPEN_MP> {
  using Base = SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::OPEN_MP>;

public:
  /**
   * @brief Forwards to SolverBase.
   * @param iters_  Number of time steps.
   * @param frames_ Number of frames to emit; @c 0 disables frame output.
   */
  OpenMPSolver(const unsigned int iters_, const unsigned int frames_)
      : Base(iters_, frames_) {};

  ~OpenMPSolver() = default;

  /**
   * @brief Runs the time loop on the host.
   *
   * Allocates the two population buffers and the velocity-norm scratch
   * array, initialises the equilibrium from the lattice's initial @c rho and
   * @c u, then iterates: fused stream-collide, buffer swap, and a frame
   * emission every @c nskips steps.
   *
   * The macroscopic fields are materialised only when they are needed -- on
   * a frame step, on the last iteration, and always on a pressure face,
   * whose boundary condition reads them. The rest of the time the moments
   * are computed and discarded, which is what keeps the pass to one read and
   * one write per population.
   *
   * @param[in,out] lattice State to evolve; @c u and @c rho are written back.
   * @param[in]     params_ Relaxation parameters.
   *
   * @note In benchmark mode (BackendProperties) no frame is emitted at all,
   *       so the measured time is the solver's and not the writer's.
   *
   * @note With @c LBM_PROFILING the per-scope timings collected by
   *       PROFILE_SCOPE are appended to the profiling CSV when the loop ends.
   */
  void solve(Lattice<dim> &lattice,
             const CollisionParams<dim, cm_t> &params_) const override {
    logging::Logger *solver_logger = logging::create_or_get_logger("solver");

    std::vector<double> ffrom(lattice.grid.getArea() * VelocitySet::ndir, 0.0);
    std::vector<double> fto(lattice.grid.getArea() * VelocitySet::ndir, 0.0);
    std::vector<float> usq(lattice.grid.getArea());

    const CollisionStrategy<dim, VelocitySet, cm_t> cs(params_);

    const auto &props =
        profiling::BackendProperties<ExecutionBackend::OPEN_MP>::get();
    const bool benchmarking = props.getBenchmarkMode();

    {
      PROFILE_SCOPE("init_equilibrium");
      init_equilibrium(ffrom, lattice);
    }

    LBM_LOG_DEBUG(solver_logger, "Equilibrium Initialized...");
    LBM_LOG_INFO(solver_logger, "System has {} logical processors.",
                 omp_get_num_procs());
    LBM_LOG_INFO(solver_logger, "The parallel section will run on {} threads.",
                 omp_get_max_threads() >= omp_get_num_procs()
                     ? omp_get_num_procs()
                     : omp_get_max_threads());

    {
      PROFILE_SCOPE("solve_total"); // wall time of the whole loop

      for (std::size_t iter = 0; iter < this->niters; iter++) {
        const bool save = !benchmarking && iter % this->nskips == 0;
        const bool store_macroscopic = save || (iter + 1 == this->niters);

        {
          PROFILE_SCOPE("stream_collide");
          update_stream_collide(ffrom, fto, usq, lattice, cs,
                                store_macroscopic);
        }
        std::swap(ffrom, fto);

        if (save) {
          write_norms(usq);
        }
      }
    }

#ifdef LBM_PROFILING
    auto &profiler = profiling::Profiler<ProfilingSchemaOpenMP>::get();
    for (const auto &[label, e] : profiling::registry()) {
      profiler.append_row(label, format::csv_format(lattice.grid.size),
                          collision_model_to_string(cm_t),
                          backend_to_string(OPEN_MP), props.getNumThreads(),
                          e.total_ms, e.total_ms / e.calls, e.calls);
    }
#endif
  }

private:

/// Equilibrium value for one direction, given local rho/u. Same formula as
/// init_equilibrium(), factored out so the sponge layer can call it per-node.
double equilibrium_i(std::size_t diridx, double r, const utils::Vector<double, dim> &u) const {
  const double u_sq = utils::ops::dot(u, u);
  const double cidotu = utils::ops::dot(VelocitySet::dir[diridx], u);
  return VelocitySet::wi[diridx] * r *
         (1.0 + numbers::invcs_2 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
}

/// 0 outside the sponge zone, rising to `max_strength` at the outflow face
/// itself. Only faces configured as OPEN_OUTFLOW absorb; walls/periodic axes
/// are left untouched (their own BC already handles them correctly).
double sponge_strength(const Grid<dim> &grid, const Solid::DomainBC<dim> &dbc,
                       const types::Coordinate<dim> &p, int width,
                       double max_strength) const {
  using utils::ops::axis;

  double s = 0.0;
  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int c = axis(p, a);

    if (dbc.low(a) == Solid::OPEN_OUTFLOW) {
      const int dist = c;                 // 0 at the face, grows inward
      if (dist < width)
        s = fmax(s, max_strength * (1.0 - static_cast<double>(dist) / width));
    }
    if (dbc.high(a) == Solid::OPEN_OUTFLOW) {
      const int dist = n - 1 - c;
      if (dist < width)
        s = fmax(s, max_strength * (1.0 - static_cast<double>(dist) / width));
    }
  }
  return s; // combine faces with max(), not sum(): avoid double-damping in corners
}

  /**
   * @brief Fills the populations with the equilibrium of the initial fields.
   *
   * @f$ f_i = f_i^{eq}(\rho_0, \mathbf{u}_0) @f$ at every node, including
   * solid ones -- their populations are never read, so leaving them at
   * equilibrium costs nothing and avoids a branch.
   *
   * @param[out] part_stream Population buffer, @c getArea()*ndir entries.
   * @param[in]  lattice     Source of the initial @c rho and @c u.
   */
  inline void init_equilibrium(std::vector<double> &part_stream,
                               const Lattice<dim> &lattice) const {
    const auto ext = lattice.grid.extents();
    const auto area = static_cast<std::ptrdiff_t>(lattice.grid.getArea());

    using utils::ops::dot;

#pragma omp parallel for shared(lattice, part_stream) schedule(static)
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

  /**
   * @brief One time step: stream, moments, equilibrium and collision, fused.
   *
   * For each fluid node, in this order:
   * 1. sum the node's own populations into @c r_wall, the density the
   *    bounce-back conditions need *before* any link is touched;
   * 2. for each direction, resolve the link (Solid::resolve_link()) and
   *    either stream from the source node or apply the boundary condition it
   *    named;
   * 3. reduce the incoming populations to @f$ \rho @f$ and @f$ \mathbf{u} @f$;
   * 4. store them if they are wanted, and the velocity norm if a frame is due;
   * 5. collide in place and write the result to @p fto.
   *
   * @param[in]     ffrom             Populations of the previous step.
   * @param[out]    fto               Populations of the new step.
   * @param[out]    usq               Velocity norms, written only when
   *                                  @p store_macroscopic is true.
   * @param[in,out] lattice           Geometry and boundary description in,
   *                                  macroscopic fields out.
   * @param[in]     cs                Collision strategy to apply.
   * @param[in]     store_macroscopic Whether @c u, @c rho and @p usq have to
   *                                  be written this step.
   *
   * @note Solid nodes are skipped on the @c solid_mask, never on a boundary
   *       type: a fluid node on a domain edge carries a face BC and must
   *       still be updated. This holds because no scheme here reads a solid
   *       node's populations; interpolated bounce-back (Bouzidi,
   *       Filippova-Hanel) or Guo extrapolation would break that assumption
   *       and require distinguishing surface from interior solid nodes.
   */
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

      // SPONGE LAYER: extra relaxation toward equilibrium near OPEN_OUTFLOW faces,
      // damping any residual wave -- acoustic or convective -- before it reaches
      // the boundary. width/max_strength are tuning knobs; start conservative.
      {
        constexpr int sponge_width = 1;
        constexpr double sponge_max = 0.3;

        const double s = sponge_strength(lattice.grid, lattice.domain_bc, p,
                                          sponge_width, sponge_max);   // <-- fix chiamata
        if (s > 0.0) {
          const double rho_ref = 1.0;
          utils::Vector<double, dim> u_ref;
          
          if constexpr (dim == 2) {
              u_ref = utils::Vector<double, 2>(0.0, 0.0);
          } else {
              u_ref = utils::Vector<double, 3>(0.0, 0.0, 0.0);
          }

          const double r_blend = r + s * (rho_ref - r);
          utils::Vector<double, dim> u_blend = u;
          u_blend.dx += s * (u_ref.dx - u.dx);
          u_blend.dy += s * (u_ref.dy - u.dy);

          for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
            fp[diridx] = equilibrium_i(diridx, r_blend, u_blend);   // <-- fix chiamata
          }
        }
      }

      // COPY LOCAL DENSITY TO GRID
#pragma omp simd
      for (std::ptrdiff_t diridx = 0;
           diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); diridx++) {
        fto[lattice.grid.field_index(p, diridx, VelocitySet::ndir)] =
            fp[diridx];
      }
    }
  };

  /**
   * @brief Emits one frame of velocity norms to the attached listeners.
   *
   * The payload is the raw image of @p usq: @c getArea() @c float32 values in
   * scalar_index() order, native byte order, with no per-frame header. See
   * the "Output formats" page.
   *
   * @param usq Velocity norms of the current step.
   */
  inline void write_norms(const std::vector<float> &usq) const {
    std::vector<char> buf(usq.size() * sizeof(float));
    std::memcpy(buf.data(), usq.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm
#endif // __LBM_SIM_SOLVER_SOLVER_2D
