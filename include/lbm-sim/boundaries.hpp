#ifndef __LBM_SIM_BOUNDARIES_HPP
#define __LBM_SIM_BOUNDARIES_HPP

#include "lbm-sim/backend.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-detection/shape.hpp"
#include "lbm-sim/constants.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/cuda/annotations.hpp"
#include "lbm-sim/types/common.hpp"

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
constexpr types::boundary_t PRESSURE_PERIODIC_INLET = 4;
constexpr types::boundary_t PRESSURE_PERIODIC_OUTLET = 5;
constexpr types::boundary_t BB_INNER = 6;

template <unsigned short int dim>
std::vector<uint8_t> compute_boundary_mask(
    const std::unordered_map<unsigned int, uint8_t> &obst_type_map,
    const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles,
    types::DimPoint<dim> size) {

  // Initialize the boundary mask with NONE for all grid points
  types::boundary_mask_t boundary_mask(utils::ops::measure(size), Solid::NONE);
  const Grid<dim> grid(size);

  // For each obstacle, get its perimeter and mark the corresponding grid points
  // in the boundary mask
  for (auto obs_idx = 0; obs_idx < obstacles.size(); obs_idx++) {
    const CollisionDetection::CollisionArea<dim> &obstacle = obstacles[obs_idx];
    const std::vector<types::Coordinate<dim>> &perimeter =
        obstacle.getPerimeter();

    // check if the obstacle is a Circle by using std::holds_alternative
    if (std::holds_alternative<lbm::CollisionDetection::Circle<dim>>(
            obstacle.collision_shapes[0])) {

      // DEBUG: stampo la class name dell'oggetto
      std::cout << "Obstacle " << obs_idx
                << " is of type: " << typeid(obstacle).name() << std::endl;

      #pragma omp parallel for shared(boundary_mask) schedule(static) collapse(2)
      // check the whole gird with contains() and mark the boundary mask for
      // each obstacle
      for (auto y = 0; y < size.y; ++y) {
        for (auto x = 0; x < size.x; ++x) {
          const types::Coordinate<dim> p(x, y);

          if (obstacle.contains(p)) {
            boundary_mask[grid.scalar_index(p)] = Solid::BB_INNER;
          }
        }
      }

      #pragma omp parallel for shared(boundary_mask, obst_type_map) schedule(static) collapse(2)
      for (auto y = 0; y < size.y; ++y) {
        for (auto x = 0; x < size.x; ++x) {

          const types::Coordinate<dim> p(x, y);
          const auto p_idx = grid.scalar_index(p);

          if (boundary_mask[p_idx] == Solid::BB_INNER) continue;

          // Now we can create the outer shell:
          for (auto diridx = 0; diridx < D2Q9::ndir; ++diridx) {
            const types::Coordinate<2> src = p - D2Q9::dir[diridx];

            // Very important: check if the src point is inside the grid, 
            // otherwise we will have an out-of-bounds access and a crash!!!
            if (!grid.contains(src)) {
              continue;
            }

            // if the point p is inside the shape we mark all its neighbors
            // which are not marked as inner as boundary nodes
            if (boundary_mask[grid.scalar_index(src)] == Solid::BB_INNER) {
              boundary_mask[p_idx] = obst_type_map.at(obs_idx);
              break; // No need to check other directions, we already marked this point
            }
          }
        }
      }

    } else {
      // Now ovverride / mark the perimeter points with the specific obstacle type
      #pragma omp parallel for shared(obst_type_map, perimeter, boundary_mask) schedule(static)
      for (auto per_idx = 0; per_idx < perimeter.size(); per_idx++) {
        const types::Coordinate<dim> &p = perimeter[per_idx];
        boundary_mask[grid.scalar_index(p)] = obst_type_map.at(obs_idx);
      }
    }
  }

  return boundary_mask;
}

template <types::dim_t dim, typename VelocitySet>
LBM_HD_FUNC inline void
apply_bb_rigid_wall(double *fp, const double *ffrom, const std::size_t diridx,
                    const Grid<dim> grid, const types::Coordinate<dim> p) {
  const auto oppdir = detail::opposite<VelocitySet>(diridx);
  fp[diridx] = ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)];
}

template <types::dim_t dim, typename VelocitySet>
LBM_HD_FUNC inline void
apply_bb_moving_wall(double *fp, const double *ffrom, const std::size_t diridx,
                     const Grid<dim> grid, const types::Coordinate<dim> p,
                     const double localrho,
                     const utils::Vector<double, dim> u0) {
  const auto oppdir = detail::opposite<VelocitySet>(diridx);
  const auto opp_direction = detail::direction<dim, VelocitySet>(oppdir);

  fp[diridx] = ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)] -
               2.0 * detail::weight<VelocitySet>(oppdir) * localrho *
                   utils::ops::dot(opp_direction, u0) * numbers::invcs_2;
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void
apply_periodic(double *fp, const double *ffrom, const std::size_t diridx,
               const Grid<dim> grid, const types::Coordinate<dim> p) {
  const auto direction = detail::direction<dim, VelocitySet>(diridx);
  const auto grid_size = static_cast<types::Coordinate<dim>>(grid.size);
  const types::Coordinate<dim> wrapped =
      (grid_size + p - direction) % grid_size;

  fp[diridx] = ffrom[grid.field_index(wrapped, diridx, VelocitySet::ndir)];
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_periodic_with_pressure_variation(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> grid, const types::Coordinate<dim> pos, const double p,
    const double *rho, const utils::Vector<double, dim> *u) {
  using utils::ops::dot;

  const auto direction = detail::direction<dim, VelocitySet>(diridx);
  const auto grid_size = static_cast<types::Coordinate<dim>>(grid.size);

  types::Coordinate<dim> wrapped = (grid_size + pos - direction) % grid_size;

  const utils::Vector<double, dim> u_wrapped = u[grid.scalar_index(wrapped)];
  const double rho_wrapped = rho[grid.scalar_index(wrapped)];

  const double omusq = -1.5 * dot(u_wrapped, u_wrapped);
  const double cidotu = dot(direction, u_wrapped);

  const auto weight = detail::weight<VelocitySet>(diridx);
  // f_i^eq(pin/pout, uwrapped) = wi * pin/pout * f(uwrapped);
  const double feq_loc =
      weight * p * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

  // fi*(wrapped,t)
  const double fi = ffrom[grid.field_index(wrapped, diridx, VelocitySet::ndir)];
  // f_i^eq(wrapped, t) = wi * rho(wrapped) * f(uwrapped);
  const double feq_from = weight * rho_wrapped *
                          (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

  fp[diridx] = feq_loc + fi - feq_from;
}

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_boundary_condition(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> grid, const types::boundary_t *boundary_mask,
    const double *rho, const utils::Vector<double, dim> *u,
    const types::Coordinate<dim> pos, const double localrho,
    const utils::Vector<double, dim> u0, const double pin, const double pout) {
  const types::boundary_t boundary = boundary_mask[grid.scalar_index(pos)];

  switch (boundary) {
  case BB_RIGID_WALL:
    apply_bb_rigid_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, pos);
    break;
  case BB_INNER:
    apply_bb_rigid_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, pos);
    break;
  case BB_MOVING_WALL:
    apply_bb_moving_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, pos,
                                           localrho, u0);
    break;
  case PERIODIC:
    apply_periodic<dim, VelocitySet>(fp, ffrom, diridx, grid, pos);
    break;

  case PRESSURE_PERIODIC_INLET:
    apply_periodic_with_pressure_variation<dim, VelocitySet>(
        fp, ffrom, diridx, grid, pos, pin, rho, u);
    break;

  case PRESSURE_PERIODIC_OUTLET:
    apply_periodic_with_pressure_variation<dim, VelocitySet>(
        fp, ffrom, diridx, grid, pos, pout, rho, u);
    break;

  default:
    break;
  }
}

} // namespace Solid
} // namespace lbm

#endif // __LBM_SIM_BOUNDARIES_HPP
