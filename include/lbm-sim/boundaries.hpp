#ifndef __LBM_SIM_BOUNDARIES_HPP
#define __LBM_SIM_BOUNDARIES_HPP

#include "collision-detection/collision-area.hpp"
#include "collision-detection/core/types.hpp"
#include "collision-detection/core/vector.hpp"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <omp.h>

namespace lbm {
namespace Solid {

namespace types {
using boundary_t = uint8_t;
using boundary_mask_t = std::vector<boundary_t>;
} // namespace types

constexpr types::boundary_t NONE = 0;
constexpr types::boundary_t BB_RIGID_WALL = 1;
constexpr types::boundary_t BB_MOVING_WALL = 2;

template <unsigned short int dim>
std::size_t coord_to_scalar(CollisionDetection::types::Coordinate<dim> p,
                            CollisionDetection::types::DimPoint<dim> size) {
  // FIXME: This is a repetition of grid.scalar_index()
  if constexpr (dim == 2) {
    return size.x * p.y + p.x;
  } else {
    static_assert(dim == 3,
                  "Solid::coord_to_scalar(): 3D Implementation Missing!");
  }
}

template <unsigned short int dim>
std::vector<uint8_t> compute_boundary_mask(
    const std::unordered_map<unsigned int, uint8_t> &strt,
    const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles,
    CollisionDetection::types::DimPoint<dim> size) {
  using CollisionDetection::types::Coordinate;

  types::boundary_mask_t boundary_mask(size.measure() * sizeof(boundary_mask),
                                       Solid::NONE);

#pragma omp parallel for shared(strt, obstacles, boundary_mask) schedule(static)
  for (std::size_t obs_idx = 0; obs_idx < obstacles.size(); obs_idx++) {
    const CollisionDetection::CollisionArea<dim> &obstacle = obstacles[obs_idx];
    for (const Coordinate<dim> &p : obstacle.getPerimeter()) {
      boundary_mask[coord_to_scalar<dim>(p, size)] = strt.at(obs_idx);
    }
  }

  return std::move(boundary_mask);
}

template <typename VelocitySet>
inline void apply_bb_rigid_wall(std::array<double, VelocitySet::ndir> &fp,
                                const unsigned int diridx) {
  fp[VelocitySet::opp[diridx]] = fp[diridx];
}

template <unsigned short int dim, typename VelocitySet>
inline void
apply_bb_moving_wall(const double &localrho,
                     const CollisionDetection::utils::Vector<double, dim> u0,
                     std::array<double, VelocitySet::ndir> &fp,
                     const unsigned int diridx) {
  // FIXME: c_s = 1/sqrt(3) should be a variable (or even
  // better 1/c_s)
  fp[VelocitySet::opp[diridx]] =
      fp[diridx] -
      2 * VelocitySet::wi[diridx] * localrho *
          CollisionDetection::utils::dot(VelocitySet::dir[diridx], u0) * 3;
}

} // namespace Solid
} // namespace lbm

#endif // __LBM_SIM_BOUNDARIES_HPP
