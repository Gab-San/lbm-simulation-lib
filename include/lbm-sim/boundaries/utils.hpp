#pragma once

#include "lbm-sim/boundaries/types.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/types/base.hpp"

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

    // popola la cache in modo sequenziale, un solo thread, nessuna race possibile
    for (const auto &shape_variant : obstacles[oid].collision_shapes) {
        std::visit([](const auto &shape) { shape.precompute(); }, shape_variant);
    }

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
  // for (types::dim_t a = 0; a < dim; ++a) {
  //   const int n = static_cast<int>(axis(grid.size, a));
  //   const int s = axis(src, a);
  //   if (s < 0 && dbc.low(a) == BB_MOVING_WALL)
  //     return {src, BB_MOVING_WALL, types::FLUID};
  //   if (s >= n && dbc.high(a) == BB_MOVING_WALL)
  //     return {src, BB_MOVING_WALL, types::FLUID};
  // }

  //     pass B: if any violated axis is an OPEN_OUTFLOW face, keep that
  //     boundary on corner/edge links before falling back to the generic
  //     x -> y -> z face selection below. This prevents a corner of an outlet
  //     from being silently classified as a wall or plain streaming.
  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int s = axis(src, a);
    if (s < 0 && dbc.low(a) == OPEN_OUTFLOW)
      return {src, OPEN_OUTFLOW, types::FLUID};
    if (s >= n && dbc.high(a) == OPEN_OUTFLOW)
      return {src, OPEN_OUTFLOW, types::FLUID};
  }

  //     pass C: first violated axis, in x -> y -> z order
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
    const auto lo = dbc.low(a);
    const auto hi = dbc.high(a);
    if (c == 0 && (is_pressure(lo) || lo == OPEN_OUTFLOW))
      return true;
    if (c == n - 1 && (is_pressure(hi) || hi == OPEN_OUTFLOW))
      return true;
  }
  return false;
}

} // namespace lbm::Solid
