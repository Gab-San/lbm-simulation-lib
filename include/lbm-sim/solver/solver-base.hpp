#pragma once

#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/data/data-observable.hpp"
#include "lbm-sim/lattice.hpp"
#include "lbm-sim/metadata.hpp"

namespace lbm {
template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t,
          enum ExecutionBackend backend_t>
class SolverBase : public DataObservable {
protected:
  const unsigned int niters, nskips;

public:
  static constexpr CollisionModel cm_type = cm_t;
  static constexpr ExecutionBackend backend_type = backend_t;

  SolverBase(const unsigned int num_iters_, const unsigned int num_frames_)
      : niters(num_iters_),
        nskips(num_frames_ > 0 ? num_iters_ / num_frames_ : num_iters_) {
    if (num_iters_ < num_frames_) {
      throw std::invalid_argument("the number of iterations must be higher "
                                  "than the number of frames!");
    }
  };

  virtual ~SolverBase() = default;
  virtual void solve(Lattice<dim> &lattice,
                     const CollisionParams<dim, cm_t> &params_) const = 0;
};
} // namespace lbm
