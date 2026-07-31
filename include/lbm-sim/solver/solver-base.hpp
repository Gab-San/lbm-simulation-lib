#pragma once

// LIB SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/backend/metadata.hpp"

#include "lbm-sim/data/data-observable.hpp"

#include "lbm-sim/structure.hpp"

// C++ STANDARD LIB
#include <vector>

namespace lbm {
template <int dim, typename VelocitySet, enum CollisionModel cm_t,
          enum ExecutionBackend backend_t>
class SolverBase : public DataObservable {
protected:
  const unsigned int niters, nskips;
  const Structure<dim> &strt;

public:
  static constexpr CollisionModel cm_type = cm_t;
  static constexpr ExecutionBackend backend_type = backend_t;

  SolverBase(const unsigned int num_iters_, const unsigned int num_frames_,
             const Structure<dim> &strt_)
      : niters(num_iters_),
        nskips(num_frames_ > 0 ? num_iters_ / num_frames_ : num_iters_),
        strt(strt_) /*, logger(logger_)*/ {
    if (num_iters_ < num_frames_) {
      throw std::invalid_argument("the number of iterations must be higher "
                                  "than the number of frames!");
    }
  };

  virtual ~SolverBase() = default;

  virtual void init_equilibrium(const Grid<dim> &grid,
                                std::vector<double> &part_stream) const = 0;
  virtual void solve(Grid<dim> &grid, const Params<dim, cm_t> &params_,
                     std::vector<double> &ffrom,
                     std::vector<double> &fto) const = 0;

protected:
  virtual void write_norms(const Grid<dim> &grid) const = 0;
};

template <enum CollisionModel cm_t, enum ExecutionBackend backend_t>
class SolverBase2D : public SolverBase<2, D2Q9, cm_t, backend_t> {
  using Base = SolverBase<2, D2Q9, cm_t, backend_t>;

public:
  SolverBase2D(const unsigned int num_iters_, const unsigned int num_frames_,
               const Structure<2> &strt_)
      : Base(num_iters_, num_frames_, strt_) {}

  virtual ~SolverBase2D() = default;
  virtual void
  init_equilibrium(const Grid<2> &grid,
                   std::vector<double> &part_stream) const override = 0;

  virtual void solve(Grid<2> &grid, const Params<2, cm_t> &params_,
                     std::vector<double> &ffrom,
                     std::vector<double> &fto) const override = 0;

protected:
  virtual void write_norms(const Grid<2> &grid) const override = 0;
};

} // namespace lbm
