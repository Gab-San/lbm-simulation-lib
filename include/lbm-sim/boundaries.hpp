#ifndef __LBM_SIM_BOUNDARIES_HPP
#define __LBM_SIM_BOUNDARIES_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/vector.hpp"

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
constexpr types::boundary_t PRESSURE_PERIODIC_INLET = 4;
constexpr types::boundary_t PRESSURE_PERIODIC_OUTLET = 5;

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

template <unsigned short int dim>
std::vector<uint8_t> compute_boundary_mask(
    const std::unordered_map<unsigned int, uint8_t> &obst_type_map,
    const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles,
    types::DimPoint<dim> size) {

  types::boundary_mask_t boundary_mask(utils::ops::measure(size), Solid::NONE);

  // FIXME: THIS SHOULD BE SUBSTITUTED BY .CONTAINS()
  // AND IT SHOULD BE NAMED COMPUTE_OBSTACLE_MASK
  for (std::size_t obs_idx = 0; obs_idx < obstacles.size(); obs_idx++) {
    const CollisionDetection::CollisionArea<dim> &obstacle = obstacles[obs_idx];
    const std::vector<types::Coordinate<dim>> &perimeter =
        obstacle.getPerimeter();
#pragma omp parallel for shared(obst_type_map, perimeter, boundary_mask)       \
    schedule(static)
    for (std::size_t per_idx = 0; per_idx < perimeter.size(); per_idx++) {
      const types::Coordinate<dim> &p = perimeter[per_idx];
      boundary_mask[coord_to_scalar<dim>(p, size)] = obst_type_map.at(obs_idx);
    }
  }

  return boundary_mask;
}

template <unsigned short int dim, typename VelocitySet>
inline void apply_bb_rigid_wall(std::array<double, VelocitySet::ndir> &fp,
                                const std::vector<double> &ffrom,
                                const std::size_t diridx, const Grid<dim> &grid,
                                const types::Coordinate<dim> p) {
  const auto oppdir = VelocitySet::opp[diridx];
  fp[diridx] = ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)];
}

template <unsigned short int dim, typename VelocitySet>
inline void
apply_bb_moving_wall(std::array<double, VelocitySet::ndir> &fp,
                     const std::vector<double> &ffrom, const std::size_t diridx,
                     const Grid<dim> &grid, const types::Coordinate<dim> p,
                     const double localrho,
                     const utils::Vector<double, dim> u0) {
  // FIXME: c_s = 1/sqrt(3) should be a variable (or even
  // better 1/c_s)
  const auto oppdir = VelocitySet::opp[diridx];
  fp[diridx] = ffrom[grid.field_index(p, oppdir, VelocitySet::ndir)] -
               2 * VelocitySet::wi[oppdir] * localrho *
                   utils::ops::dot(VelocitySet::dir[oppdir], u0) * 3;
}

template <unsigned short int dim, typename VelocitySet>
inline void apply_periodic(std::array<double, VelocitySet::ndir> &fp,
                           const std::vector<double> &ffrom,
                           const std::size_t diridx, const Grid<dim> &grid,
                           const types::Coordinate<dim> p) {
  types::Coordinate<dim> wrapped =
      (static_cast<types::Coordinate<dim>>(grid.size) + p -
       VelocitySet::dir[diridx]) %
      static_cast<types::Coordinate<dim>>(grid.size);
  fp[diridx] = ffrom[grid.field_index(wrapped, diridx, VelocitySet::ndir)];
}

template <unsigned short int dim, typename VelocitySet>
inline void apply_periodic_with_pressure_variation(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> grid, const types::Coordinate<dim> pos, const double *rho,
    const double p, const utils::Vector<double, dim> *u) {
  using utils::ops::dot;

  types::Coordinate<dim> wrapped =
      (static_cast<types::Coordinate<dim>>(grid.size) + pos -
       VelocitySet::dir[diridx]) %
      static_cast<types::Coordinate<dim>>(grid.size);

  const utils::Vector<double, dim> u_wrapped = u[grid.scalar_index(wrapped)];
  const double rho_wrapped = rho[grid.scalar_index(wrapped)];

  const double omusq = -1.5 * dot(u_wrapped, u_wrapped);
  const double cidotu = dot(VelocitySet::dir[diridx], u_wrapped);

  // f_i^eq(pin/pout, uwrapped) = wi * pin/pout * f(uwrapped);
  const double feq_loc = VelocitySet::wi[diridx] * p *
                         (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

  // fi*(wrapped,t)
  const double fi = ffrom[grid.field_index(wrapped, diridx, VelocitySet::ndir)];
  // f_i^eq(wrapped, t) = wi * rho(wrapped) * f(uwrapped);
  const double feq_from = VelocitySet::wi[diridx] * rho_wrapped *
                          (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

  fp[diridx] = feq_loc + fi - feq_from;
}

} // namespace Solid
} // namespace lbm

#endif // __LBM_SIM_BOUNDARIES_HPP
