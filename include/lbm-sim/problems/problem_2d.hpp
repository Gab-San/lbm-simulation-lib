#pragma once

#include "lbm-sim/lattice.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

#include <omp.h>

namespace lbm {

class LidCavity2D {
public:
  virtual ~LidCavity2D() = default;

  void
  init(Lattice<2> &lattice, const utils::Vector<double, 2> &init_vel,
       const std::vector<types::Coordinate<2>> &initialization_points) const {
    // the problem shall be set up by the CPU
    // set lid velocity at the top boundary
    for (const auto &p : initialization_points) {
      lattice.u[lattice.grid.scalar_index(p)] = init_vel;
    }
  }
};

} // namespace lbm
