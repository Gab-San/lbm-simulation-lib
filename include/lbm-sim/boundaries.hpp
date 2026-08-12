#ifndef __LBM_SIM_BOUNDARIES_HPP
#define __LBM_SIM_BOUNDARIES_HPP

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/collision-detection/collision-area.hpp"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <omp.h>

namespace lbm {

namespace types {
using boundary_t = uint8_t;
using boundary_mask_t = std::vector<boundary_t>;
} // namespace types

namespace Solid {

constexpr types::boundary_t NONE = 0;
constexpr types::boundary_t BB_RIGID_WALL = 1;
constexpr types::boundary_t BB_MOVING_WALL = 2;
constexpr types::boundary_t PERIODIC = 3;

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

template <typename VelocitySet>
inline void apply_bb_rigid_wall(std::array<double, VelocitySet::ndir> &fp,
                                const unsigned int diridx) {
  fp[VelocitySet::opp[diridx]] = fp[diridx];
}

template <unsigned short int dim, typename VelocitySet>
inline void apply_bb_moving_wall(const double &localrho,
                                 const utils::Vector<double, dim> u0,
                                 std::array<double, VelocitySet::ndir> &fp,
                                 const unsigned int diridx) {
  // FIXME: c_s = 1/sqrt(3) should be a variable (or even
  // better 1/c_s)
  fp[VelocitySet::opp[diridx]] =
      fp[diridx] - 2 * VelocitySet::wi[diridx] * localrho *
                       utils::ops::dot(VelocitySet::dir[diridx], u0) * 3;
}

// FIXME: Make dimension independent

// TODO: Check implementation
template <unsigned short int dim, typename VelocitySet>
inline void apply_continue(const double &localrho, const Grid<2> &grid,
                           const std::vector<double> &ffrom,
                           types::Coordinate<2> p,
                           std::array<double, VelocitySet::ndir> &fp,
                           const unsigned int diridx) {
  using Coord = types::Coordinate<2>;

  const unsigned int missing = VelocitySet::opp[diridx];

  const Coord raw = p + VelocitySet::dir[diridx];
  const int size_x = static_cast<int>(grid.size.x);
  const int wrapped_x = ((raw.x % size_x) + size_x) % size_x;
  const Coord source(wrapped_x, raw.y);

  fp[missing] = ffrom[grid.field_index(source, missing, VelocitySet::ndir)];
}

// NOTE: Why is this all commented?
/*const types::Coordinate<2> x=p -VelocitySet::dir[diridx];
using Coord = types::Coordinate<2>;

if(p.x==grid.size.x-1){
  if (p.y == grid.size.y - 1) {
    if(diridx==5)fp[8] = fp[diridx];
  }else if (p.y == 0) {
    if(diridx==8)fp[5] = fp[diridx];
  } else {
    // left boundary

    if(diridx==1)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x+Coord(grid.size.x-1,0),
diridx,VelocitySet::ndir )]; if(diridx==5)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x+Coord(grid.size.x-1,1),
diridx,VelocitySet::ndir )]; if(diridx==8)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x+Coord(grid.size.x-1,-1),
diridx,VelocitySet::ndir )];

  }
}else{
  if (p.y == grid.size.y - 1) {
    if(diridx==6)fp[7] = fp[diridx];
  }else if (p.y == 0) {
    if(diridx==7)fp[6] = fp[diridx];
  } else {
    // left boundary
    if(diridx==3)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x-Coord(grid.size.x-1,0),
diridx,VelocitySet::ndir )]; if(diridx==6)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x-Coord(grid.size.x-1,1),
diridx,VelocitySet::ndir )]; if(diridx==7)ffrom[grid.field_index(x, diridx,
VelocitySet::ndir)]=ffrom[grid.field_index(x-Coord(grid.size.x-1,-1),
diridx,VelocitySet::ndir )];
  }
}

*/

} // namespace Solid
} // namespace lbm

#endif // __LBM_SIM_BOUNDARIES_HPP
