#pragma once

#include "lbm-sim/backend.hpp"
#include "lbm-sim/boundaries/types.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/constants.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/cuda/annotations.hpp"
#include "lbm-sim/types/common.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

#include <omp.h>

namespace lbm::Solid {

/**
 * Paint one obstacle id per solid node. Runs once, at setup.
 *
 * Unlike the old compute_boundary_mask() this stores an **id**, not a BC type,
 * and loops over the clamped AABB of each shape rather than over the grid.
 * Storing the id is what lets a per-obstacle wall velocity reach
 * apply_bb_moving_wall().
 *
 * Race-free: schedule(static) partitions the outermost index, so two threads
 * never touch the same row. Overlapping obstacles are last-writer-wins across
 * the sequential `oid` loop, which is deterministic. This runs once at setup,
 * so the per-node contains() cost is irrelevant.
 */
template <types::dim_t dim>
types::solid_mask_t compute_solid_mask(
    const std::vector<CollisionDetection::CollisionArea<dim>> &obstacles,
    const types::DimPoint<dim> size) {

  using utils::ops::axis;

  assert(obstacles.size() < types::FLUID && "obstacle id overflow");

  types::solid_mask_t mask(utils::ops::measure(size), types::FLUID);
  const Grid<dim> grid(size);

  for (std::size_t oid = 0; oid < obstacles.size(); ++oid) {
    const auto box = obstacles[oid].aabb();

    types::Coordinate<dim> lo = box.min;
    types::Coordinate<dim> hi = box.max;

    bool empty = false;
    for (types::dim_t a = 0; a < dim; ++a) {
      const int n = static_cast<int>(axis(size, a));
      int &l = axis(lo, a);
      int &h = axis(hi, a);
      if (l < 0)
        l = 0;
      if (h > n - 1)
        h = n - 1;
      if (l > h)
        empty = true; // shape entirely outside the domain
    }
    if (empty)
      continue;

    const auto id = static_cast<types::obstacle_id_t>(oid);

    if constexpr (dim == 2) {
#pragma omp parallel for shared(mask, obstacles, grid) schedule(static)
      for (int y = lo.y; y <= hi.y; ++y) {
        for (int x = lo.x; x <= hi.x; ++x) {
          const types::Coordinate<2> q(x, y);
          if (obstacles[oid].contains(q))
            mask[grid.scalar_index(q)] = id;
        }
      }
    } else {
#pragma omp parallel for shared(mask, obstacles, grid) schedule(static)
      for (int z = lo.z; z <= hi.z; ++z) {
        for (int y = lo.y; y <= hi.y; ++y) {
          for (int x = lo.x; x <= hi.x; ++x) {
            const types::Coordinate<3> q(x, y, z);
            if (obstacles[oid].contains(q))
              mask[grid.scalar_index(q)] = id;
          }
        }
      }
    }
  }
  return mask;
}

/// What a single link (p, dir) resolves to.
template <types::dim_t dim> struct LinkResolution {
  types::Coordinate<dim> src; // already wrapped where applicable
  types::boundary_t bc;       // NONE -> plain streaming from src
  types::obstacle_id_t oid;   // FLUID unless a solid owns the link
};

/**
 * The whole boundary decision for one link, in one function: replaces the
 * `grid.contains(src)` test in the solver, apply_periodic(), and the mask
 * lookup that used to live inside apply_boundary_condition().
 *
 * A single mask slot per node cannot express this. Take a 1000^2 grid,
 * periodic in x, rigid bottom, p = (0,0): with dir = (1,1) the source
 * (-1,-1) wraps in x to (999,-1), then step 2 sees y < 0 and hands the link
 * to the bottom wall; with dir = (1,0) it wraps to (999,0), which is fluid,
 * and the link streams normally. Two directions, one node, two answers.
 *
 * PRIORITY CHOICE at step 2, recorded because the Ghia comparison is
 * sensitive to it at the two top corners of a lid-driven cavity: when a
 * diagonal link crosses two non-periodic faces, a **moving wall wins** over a
 * rigid one, whichever axis it sits on. This reproduces the pre-refactor
 * behaviour, where the lid was painted into the boundary mask after the side
 * walls and so owned the corner nodes. Everything else falls back to axis
 * order x -> y -> z.
 */
template <types::dim_t dim>
LBM_HD_FUNC inline LinkResolution<dim> resolve_link(
    const Grid<dim> &grid, const DomainBC<dim> &dbc,
    const types::obstacle_id_t *solid_mask, const ObstacleData<dim> *obstacles,
    const types::Coordinate<dim> &p, const types::Direction<dim> &dir) {

  using utils::ops::axis;

  types::Coordinate<dim> src = p - dir;
  types::boundary_t wrap_bc = NONE;

  // (1) wrap every axis whose face is periodic.
  //     |dir| <= 1 per component, so one add/sub is always enough.
  for (types::dim_t a = 0; a < dim; ++a) {
    if (!dbc.is_periodic_axis(a))
      continue;
    const int n = static_cast<int>(axis(grid.size, a));
    int &s = axis(src, a);
    if (s < 0) {
      s += n;
      if (is_pressure(dbc.low(a)))
        wrap_bc = dbc.low(a);
    } else if (s >= n) {
      s -= n;
      if (is_pressure(dbc.high(a)))
        wrap_bc = dbc.high(a);
    }
  }

  // (2) still outside? the responsible non-periodic face owns the link.
  //     A wall beats a wrap, because this runs after step 1.
  //
  //     pass A: a moving wall claims the link if any violated face carries one
  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int s = axis(src, a);
    if (s < 0 && dbc.low(a) == BB_MOVING_WALL)
      return {src, BB_MOVING_WALL, types::FLUID};
    if (s >= n && dbc.high(a) == BB_MOVING_WALL)
      return {src, BB_MOVING_WALL, types::FLUID};
  }

  //     pass B: first violated axis, in x -> y -> z order
  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int s = axis(src, a);
    if (s < 0)
      return {src, dbc.low(a), types::FLUID};
    if (s >= n)
      return {src, dbc.high(a), types::FLUID};
  }

  // (3) inside the domain: is the source node solid?
  const types::obstacle_id_t oid = solid_mask[grid.scalar_index(src)];
  if (oid != types::FLUID)
    return {src, obstacles[oid].bc_type, oid};

  // (4) fluid. wrap_bc is NONE for ordinary streaming, or the pressure face
  //     whose rescale still has to be applied to the already-wrapped src.
  return {src, wrap_bc, types::FLUID};
}

/// Pure coordinate test replacing the macroscopic-store mask read in the
/// solver: is `p` itself sitting on a pressure-driven domain face?
template <types::dim_t dim>
LBM_HD_FUNC inline bool on_pressure_face(const Grid<dim> &grid,
                                         const DomainBC<dim> &dbc,
                                         const types::Coordinate<dim> &p) {
  using utils::ops::axis;
  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int c = axis(p, a);
    if (c == 0 && is_pressure(dbc.low(a)))
      return true;
    if (c == n - 1 && is_pressure(dbc.high(a)))
      return true;
  }
  return false;
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

/// Pressure-periodic rescale on an already-wrapped source node. resolve_link()
/// has done the wrap, so nothing is recomputed here.
template <types::dim_t dim, typename VelocitySet>
LBM_HD_FUNC inline void
apply_pressure_rescale(double *fp, const double *ffrom,
                       const std::size_t diridx, const Grid<dim> grid,
                       const types::Coordinate<dim> src, const double p_target,
                       const double *rho, const utils::Vector<double, dim> *u) {
  using utils::ops::dot;

  const auto direction = detail::direction<dim, VelocitySet>(diridx);
  const std::size_t si = grid.scalar_index(src); // computed once

  const utils::Vector<double, dim> u_src = u[si];
  const double rho_src = rho[si];

  const double omusq = -1.5 * dot(u_src, u_src);
  const double cidotu = dot(direction, u_src);
  const double weight = detail::weight<VelocitySet>(diridx);
  const double shape = 1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq;

  // f_i^eq(p_target, u_src) and f_i^eq(rho_src, u_src)
  const double feq_loc = weight * p_target * shape;
  const double feq_from = weight * rho_src * shape;
  const double fi = ffrom[grid.field_index(src, diridx, VelocitySet::ndir)];

  fp[diridx] = feq_loc + fi - feq_from;
}

/// Dispatch on an already-resolved link: no mask read, no re-derivation of src.
template <types::dim_t dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_boundary_condition(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> grid, const LinkResolution<dim> &link,
    const ObstacleData<dim> *obstacles, const double *rho,
    const utils::Vector<double, dim> *u, const types::Coordinate<dim> pos,
    const double localrho, const utils::Vector<double, dim> u_domain,
    const double pin, const double pout) {

  switch (link.bc) {
  case BB_RIGID_WALL:
    apply_bb_rigid_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, pos);
    break;

  case BB_MOVING_WALL: {
    // Per-obstacle wall velocity, falling back to the domain value for a face
    // BC. Storing an id rather than a BC type in the mask is precisely what
    // makes this branch possible: before, apply_bb_moving_wall() could only
    // ever receive the single global cs.params.init_vel.
    const utils::Vector<double, dim> u0 =
        (link.oid != types::FLUID) ? obstacles[link.oid].wall_velocity
                                   : u_domain;
    apply_bb_moving_wall<dim, VelocitySet>(fp, ffrom, diridx, grid, pos,
                                           localrho, u0);
    break;
  }

  case PRESSURE_PERIODIC_INLET:
    apply_pressure_rescale<dim, VelocitySet>(fp, ffrom, diridx, grid, link.src,
                                             pin, rho, u);
    break;

  case PRESSURE_PERIODIC_OUTLET:
    apply_pressure_rescale<dim, VelocitySet>(fp, ffrom, diridx, grid, link.src,
                                             pout, rho, u);
    break;

  default:
    break; // NONE never reaches here -- the caller streams instead
  }
}

} // namespace lbm::Solid
