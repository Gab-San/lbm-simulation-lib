#ifndef __LBM_SIM_STRUCTURE_HPP
#define __LBM_SIM_STRUCTURE_HPP

// COLLISION DETECTION LIB
#include "collision-detection/collision-area.hpp"

template <int dim> struct Structure {
  const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles;
  const std::size_t moving_boundary;
  Structure(const std::vector<CollisionDetection::CollisionArea<2>> &obstacles_,
            const std::size_t moving_boundary_)
      : obstacles(obstacles_), moving_boundary(moving_boundary_) {}
};

#endif // __LBM_SIM_STRUCTURE_HPP
