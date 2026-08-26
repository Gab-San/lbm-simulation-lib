#pragma once

#include "lbm-sim/types/base.hpp"
#include <vector>

namespace lbm {
namespace types {

using solid_mask_t = std::vector<obstacle_id_t>;

/// Identifies a fluid node.
inline constexpr obstacle_id_t FLUID = 0xFFFF;

} // namespace types

namespace Solid {

constexpr types::boundary_t NONE = 0;
constexpr types::boundary_t BB_RIGID_WALL = 1;
constexpr types::boundary_t BB_MOVING_WALL = 2;
constexpr types::boundary_t PERIODIC = 3;
constexpr types::boundary_t PRESSURE_PERIODIC_INLET = 4;
constexpr types::boundary_t PRESSURE_PERIODIC_OUTLET = 5;

} // namespace Solid

} // namespace lbm
