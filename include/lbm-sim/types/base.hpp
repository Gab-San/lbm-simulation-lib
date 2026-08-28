#pragma once

#include <cstdint>

namespace lbm::types {

using dim_t = unsigned short int;

using boundary_t = uint8_t;

/// Index into Lattice::obstacles. types::FLUID marks a fluid node.
/// See boundaries/types.hpp for the sentinel and for solid_mask_t.
using obstacle_id_t = uint16_t;

} // namespace lbm::types
