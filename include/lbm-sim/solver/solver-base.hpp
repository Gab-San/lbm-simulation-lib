/**
 * @file solver-base.hpp
 * @brief SolverBase: the interface LBMSimulation::solve() accepts, and the
 *        frame bookkeeping every backend shares.
 *
 * The base owns two things and no more: the iteration count and the frame
 * stride derived from it. Populations, their layout and the time loop belong
 * to the concrete backend.
 *
 * It is also a DataObservable, and a *different* one from the simulation:
 * the simulation emits the grid header, the solver emits the frames. A
 * listener that needs both has to be attached to both.
 */

#pragma once

#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/data/data-observable.hpp"
#include "lbm-sim/lattice.hpp"
#include "lbm-sim/metadata.hpp"

namespace lbm {

/**
 * @brief Common base of the solver backends.
 *
 * @tparam dim         Spatial dimension (2 or 3).
 * @tparam VelocitySet Discrete velocity set.
 * @tparam cm_t        Collision model.
 * @tparam backend_t   Execution backend. It is part of the type so that
 *                     LBMSimulation::solve() can deduce it, which is what
 *                     lets one simulation type accept either backend.
 */
template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t,
          enum ExecutionBackend backend_t>
class SolverBase : public DataObservable {
protected:
  /// Total number of time steps, and the stride between two emitted frames.
  const unsigned int niters, nskips;

public:
  /// The collision model, readable from generic code.
  static constexpr CollisionModel cm_type = cm_t;

  /// The execution backend, readable from generic code.
  static constexpr ExecutionBackend backend_type = backend_t;

  /**
   * @brief Fixes the iteration count and derives the frame stride.
   *
   * @param num_iters_  Number of time steps to run.
   * @param num_frames_ Number of frames to emit over the run. @c 0 means no
   *                    frames, in which case @c nskips is set to
   *                    @p num_iters_ so the modulo in the loop never fires
   *                    after step 0.
   *
   * @throws std::invalid_argument if @p num_frames_ exceeds @p num_iters_,
   *         which would ask for more frames than there are steps.
   *
   * @note The stride is integer division, so @c num_iters_ % @c num_frames_
   *       steps at the tail may emit one frame more or less than requested.
   */
  SolverBase(const unsigned int num_iters_, const unsigned int num_frames_)
      : niters(num_iters_),
        nskips(num_frames_ > 0 ? num_iters_ / num_frames_ : num_iters_) {
    if (num_iters_ < num_frames_) {
      throw std::invalid_argument("the number of iterations must be higher "
                                  "than the number of frames!");
    }
  };

  virtual ~SolverBase() = default;

  /**
   * @brief Runs the time loop.
   *
   * @param[in,out] lattice State to evolve. The solver reads the geometry,
   *                        the mask and the boundary conditions, and writes
   *                        the macroscopic fields @c u and @c rho back at the
   *                        end of the run (and on every frame step).
   * @param[in]     params_ Relaxation parameters for @p cm_t.
   *
   * @c const because the solver holds no mutable state of its own: the
   * populations are local to the call, so the same solver object can be run
   * on several lattices in turn.
   */
  virtual void solve(Lattice<dim> &lattice,
                     const CollisionParams<dim, cm_t> &params_) const = 0;
};
} // namespace lbm
