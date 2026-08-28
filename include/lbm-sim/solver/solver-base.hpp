#ifndef __LBM_SIM_SOLVER_SOLVER_BASE_HPP
#define __LBM_SIM_SOLVER_SOLVER_BASE_HPP

// LIB SIM LIB
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
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

template <enum CollisionModel cm_t, enum ExecutionBackend backend_t>
class SolverBase2D : public SolverBase<2, D2Q9, cm_t, backend_t> {
  using Base = SolverBase<2, D2Q9, cm_t, backend_t>;

public:
  SolverBase2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {}

  virtual ~SolverBase2D() = default;

  virtual void
  solve(Lattice<2> &lattice,
        const CollisionParams<2, cm_t> &params_) const override = 0;
};

template <typename VelocitySet, enum CollisionModel cm_t,
          enum ExecutionBackend backend_t>
class SolverBase3D : public SolverBase<3, VelocitySet, cm_t, backend_t> {
  static_assert(VelocitySet::dim == 3,
                "SolverBase3D richiede un VelocitySet con dim == 3 "
                "(es. D3Q27, D3Q19)");
  using Base = SolverBase<3, VelocitySet, cm_t, backend_t>;

public:
  SolverBase3D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {}

  virtual ~SolverBase3D() = default;

  virtual void
  solve(Lattice<3> &lattice,
        const CollisionParams<3, cm_t> &params_) const override = 0;
};

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_BASE_HPP
