/**
 * @file boundary-conditions.hpp
 * @brief The boundary kernels themselves, applied to an already-resolved
 *        link.
 *
 * By the time these run, resolve_link() has already wrapped the coordinate,
 * identified the owning obstacle and named the condition; nothing here reads
 * the mask or re-derives a source node. The division of labour is
 * deliberate: one decision per link in boundaries/utils.hpp, one arithmetic
 * kernel per condition here.
 *
 * Implemented conditions:
 * - halfway bounce-back off a stationary wall, @f$ f_i \leftarrow
 *   f_{\bar\imath} @f$;
 * - halfway bounce-back off a moving wall, with the standard momentum
 *   correction @f$ -2 w_{\bar\imath} \rho\, (\mathbf{c}_{\bar\imath} \cdot
 *   \mathbf{u}_w) / c_s^2 @f$;
 * - pressure-periodic rescale, which streams from the wrapped source but
 *   substitutes the target pressure into the equilibrium part.
 *
 * All kernels are @c LBM_HD_FUNC and shared by both backends.
 */

#pragma once

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/backend/utils.hpp"
#include "lbm-sim/boundaries/types.hpp"
#include "lbm-sim/constants.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/types/common.hpp"

#include <cassert>
#include <cstddef>

namespace lbm::Solid {

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

template <types::dim_t dim, typename VelocitySet>
LBM_HD_FUNC inline void apply_open_outflow(
    double *fp, const double *ffrom, const std::size_t diridx,
    const Grid<dim> &grid, const types::Coordinate<dim> &p,
    const types::Coordinate<dim> &src) {

  using utils::ops::axis;

  types::Coordinate<dim> upstream = p;

  for (types::dim_t a = 0; a < dim; ++a) {
    const int n = static_cast<int>(axis(grid.size, a));
    const int s = axis(src, a);
    if (s < 0) {
      axis(upstream, a) += 1; // low face: step toward +axis
      break;
    }
    if (s >= n) {
      axis(upstream, a) -= 1; // high face: step toward -axis
      break;
    }
  }

  fp[diridx] = ffrom[grid.field_index(upstream, diridx, VelocitySet::ndir)];
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

  case OPEN_OUTFLOW:
    apply_open_outflow<dim, VelocitySet>(fp, ffrom, diridx, grid, pos, link.src);
    break;

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
