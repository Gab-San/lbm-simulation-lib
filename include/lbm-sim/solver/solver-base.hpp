#pragma once

// LIB SIM LIB
#include "lbm-sim/backend/metadata.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/data-observable.hpp"

#include "lbm-sim/lattice.hpp"

// C++ STANDARD LIB
#include <vector>

namespace lbm {
template <unsigned short int dim, typename VelocitySet,
          enum CollisionModel cm_t, enum ExecutionBackend backend_t>
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

  virtual void init_equilibrium(const Lattice<dim> &lattice,
                                std::vector<double> &part_stream) const = 0;
  virtual void solve(Lattice<dim> &lattice, const Params<dim, cm_t> &params_,
                     std::vector<double> &ffrom,
                     std::vector<double> &fto) const = 0;

protected:
  virtual void write_norms(const Lattice<dim> &lattice) const = 0;
};

template <enum CollisionModel cm_t, enum ExecutionBackend backend_t>
class SolverBase2D : public SolverBase<2, D2Q9, cm_t, backend_t> {
  using Base = SolverBase<2, D2Q9, cm_t, backend_t>;

public:
  SolverBase2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {}

  virtual ~SolverBase2D() = default;
  virtual void
  init_equilibrium(const Lattice<2> &lattice,
                   std::vector<double> &part_stream) const override = 0;

  virtual void solve(Lattice<2> &lattice, const Params<2, cm_t> &params_,
                     std::vector<double> &ffrom,
                     std::vector<double> &fto) const override = 0;

protected:
  virtual void write_norms(const Lattice<2> &grid) const override = 0;
};

} // namespace lbm
