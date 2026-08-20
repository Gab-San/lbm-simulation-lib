#ifndef __LBM_SIM_BOUNDARIES_HPP
#define __LBM_SIM_BOUNDARIES_HPP

#include "lbm-sim/backend/cuda-annotations.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/collision-detection/collision-area.hpp"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <omp.h>

namespace lbm {

namespace Solid {

constexpr types::boundary_t NONE = 0;
constexpr types::boundary_t BB_RIGID_WALL = 1;
constexpr types::boundary_t BB_MOVING_WALL = 2;
constexpr types::boundary_t PERIODIC = 3;
//constexpr types::boundary_t PRESSURE_PERIODIC_INLET = 4;
//constexpr types::boundary_t PRESSURE_PERIODIC_OUTLET = 5;

template <unsigned short int dim>
std::size_t coord_to_scalar(types::Coordinate<dim> p,
                            types::DimPoint<dim> size) {
  // FIXME: This is a repetition of grid.scalar_index()
  if constexpr (dim == 2) {
    return size.x * p.y + p.x;
  } else {
    static_assert(assertion::always_false<dim>,
                  "Solid::coord_to_scalar(): 3D Implementation Missing!");
  }
}

namespace detail { //using cuda or openmp

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline types::VectorInt<dim> direction(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_dir<dim, VelocitySet>[diridx];
#else
  return VelocitySet::dir[diridx];
#endif
}

template <typename VelocitySet>
LBM_HD_FUNC inline double weight(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_wi<VelocitySet>[diridx];
#else
  return VelocitySet::wi[diridx];
#endif
}

template <typename VelocitySet>
LBM_HD_FUNC inline std::size_t opposite(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_opp<VelocitySet>[diridx];
#else
  return VelocitySet::opp[diridx];
#endif
}

} // namespace detail


template <unsigned short int dim>
std::vector<uint8_t> compute_boundary_mask(
    const std::unordered_map<unsigned int, uint8_t> &obst_type_map,
    const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles,
    types::DimPoint<dim> size) {

  types::boundary_mask_t boundary_mask(utils::ops::measure(size), Solid::NONE);
  const Grid<dim> grid(size);

  // FIXME: THIS SHOULD BE SUBSTITUTED BY .CONTAINS()
  // AND IT SHOULD BE NAMED COMPUTE_OBSTACLE_MASK
  //use different approach, matching streaming support for internal solid nodes
  for (std::size_t obs_idx = 0; obs_idx < obstacles.size(); obs_idx++) {
    const CollisionDetection::CollisionArea<dim> &obstacle = obstacles[obs_idx];
    const std::vector<types::Coordinate<dim>> &perimeter =
        obstacle.getPerimeter();
#pragma omp parallel for shared(obst_type_map, perimeter, boundary_mask)       \
    schedule(static)
    for (std::size_t per_idx = 0; per_idx < perimeter.size(); per_idx++) {
      const types::Coordinate<dim> &p = perimeter[per_idx];
      boundary_mask[grid.scalar_index(p)] = obst_type_map.at(obs_idx);
    }
  }

  return boundary_mask;
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_bb_rigid_wall(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> &grid, const types::Coordinate<dim> p) {
  const auto oppdir = detail::opposite<VelocitySet>(diridx);
  fp[diridx] = ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)];
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_bb_moving_wall(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> &grid, const types::Coordinate<dim> p,
    const double localrho, const utils::Vector<double, dim> u0) {
  const auto oppdir = detail::opposite<VelocitySet>(diridx);
  const auto opp_direction = detail::direction<dim, VelocitySet>(oppdir);

  fp[diridx] =
      ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)] -
      2.0 * detail::weight<VelocitySet>(oppdir) * localrho *
          utils::ops::dot(opp_direction, u0) * VelocitySet::inv_cs2;
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_periodic(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> &grid, const types::Coordinate<dim> p) {
  const auto direction = detail::direction<dim, VelocitySet>(diridx);
  const auto grid_size = static_cast<types::Coordinate<dim>>(grid.size);
  const types::Coordinate<dim> wrapped = (grid_size + p - direction) % grid_size;

  fp[diridx] = ffrom[grid.field_index(wrapped, diridx, VelocitySet::ndir)];
}

//cuda or openmp implem 
template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_boundary_condition(
    const types::boundary_t *boundary_mask, double *fp, const double *ffrom,
    const std::size_t diridx, const Grid<dim> &grid,
    const types::Coordinate<dim> p, const double localrho,
    const utils::Vector<double, dim> u0) {
  const types::boundary_t boundary = boundary_mask[grid.scalar_index(p)];

  switch (boundary) {
  case Solid::BB_RIGID_WALL:
    apply_bb_rigid_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, p);
    break;
  case Solid::BB_MOVING_WALL:
    apply_bb_moving_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, p, localrho,
                                           u0);
    break;
  case Solid::PERIODIC:
    apply_periodic<dim, VelocitySet>(fp, ffrom, diridx, grid, p);
    break;
  default:
    break;
  }
}

} // namespace Solid
} // namespace lbm

#endif // __LBM_SIM_BOUNDARIES_HPP
