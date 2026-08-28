#pragma once

#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/cuda/annotations.hpp"
#include "lbm-sim/types/base.hpp"
#include "lbm-sim/types/common.hpp"

#include <cassert>
#include <vector>

namespace lbm {
namespace types {

/// One id per node: which obstacle owns it, or FLUID.
///
/// Two bytes per node is 1.9 MiB on a 1000^2 grid against ~137 MiB of
/// populations, i.e. 1.4%. Not a constraint.
using solid_mask_t = std::vector<obstacle_id_t>;

/// Identifies a fluid node.
inline constexpr obstacle_id_t FLUID = 0xFFFF;

} // namespace types

namespace Solid {

constexpr types::boundary_t NONE = 0; // plain streaming
constexpr types::boundary_t BB_RIGID_WALL = 1;
constexpr types::boundary_t BB_MOVING_WALL = 2;
constexpr types::boundary_t PERIODIC = 3;
constexpr types::boundary_t PRESSURE_PERIODIC_INLET = 4;
constexpr types::boundary_t PRESSURE_PERIODIC_OUTLET = 5;

LBM_HD_FUNC constexpr bool is_pressure(const types::boundary_t b) {
  return b == PRESSURE_PERIODIC_INLET || b == PRESSURE_PERIODIC_OUTLET;
}

/// Any BC whose link resolution starts with a coordinate wrap.
LBM_HD_FUNC constexpr bool wraps(const types::boundary_t b) {
  return b == PERIODIC || is_pressure(b);
}

/**
 * Boundary condition of the six (four in 2D) domain faces.
 * O(1) in the grid size: 4 bytes in 2D, 6 in 3D, for any resolution.
 * Trivially copyable -> passed by value as a CUDA kernel argument.
 *
 * A raw C array, not std::array: std::array::operator[] is not
 * __host__ __device__ unless built with --expt-relaxed-constexpr, and a raw
 * array is equally trivially copyable.
 */
template <types::dim_t dim> struct DomainBC {
  // layout: [low(0), high(0), low(1), high(1), ...]
  types::boundary_t face[2 * dim] = {};

  LBM_HD_FUNC constexpr types::boundary_t &low(types::dim_t a) {
    return face[2 * a];
  }
  LBM_HD_FUNC constexpr types::boundary_t &high(types::dim_t a) {
    return face[2 * a + 1];
  }
  LBM_HD_FUNC constexpr types::boundary_t low(types::dim_t a) const {
    return face[2 * a];
  }
  LBM_HD_FUNC constexpr types::boundary_t high(types::dim_t a) const {
    return face[2 * a + 1];
  }

  /// True when links crossing this axis must be wrapped before anything else.
  /// Note it is true if *either* face wraps -- see
  /// assert_consistent_domain_bc().
  LBM_HD_FUNC constexpr bool is_periodic_axis(types::dim_t a) const {
    return wraps(low(a)) || wraps(high(a));
  }
};

/// Setup-time sanity check: a periodic axis that wraps on one face and walls
/// on the other is legal here but almost never what was meant.
template <types::dim_t dim>
inline void assert_consistent_domain_bc(const DomainBC<dim> &dbc) {
  for (types::dim_t a = 0; a < dim; ++a) {
    assert(wraps(dbc.low(a)) == wraps(dbc.high(a)) &&
           "periodic axis must wrap on both faces");
  }
  (void)dbc;
}

/// Side table: obstacle id -> what that obstacle does to a link.
template <types::dim_t dim> struct ObstacleData {
  types::boundary_t bc_type = BB_RIGID_WALL;
  utils::Vector<double, dim> wall_velocity{};
};

/// What a single link (p, dir) resolves to.
template <types::dim_t dim> struct LinkResolution {
  types::Coordinate<dim> src; // already wrapped where applicable
  types::boundary_t bc;       // NONE -> plain streaming from src
  types::obstacle_id_t oid;   // FLUID unless a solid owns the link
};

} // namespace Solid

} // namespace lbm
