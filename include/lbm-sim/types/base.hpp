/**
 * @file base.hpp
 * @brief The scalar type aliases the whole library is parameterised on.
 *
 * Deliberately dependency-free (only @c <cstdint>), so it can be included
 * from anywhere -- including headers that must stay clear of @c <omp.h> and
 * of the CUDA runtime -- without dragging anything else in.
 *
 * The widths are chosen for what they cost per node, not for convenience:
 * see boundaries/types.hpp for the arithmetic on the mask.
 */

#pragma once

#include <cstdint>

namespace lbm::types {

/// Spatial dimension, 2 or 3. A non-type template parameter almost
/// everywhere, hence the small unsigned type.
using dim_t = unsigned short int;

/// A boundary condition tag; the values are the constants in
/// @c lbm::Solid (NONE, BB_RIGID_WALL, ...). One byte, because a DomainBC
/// holds @c 2*dim of them and is passed by value to a CUDA kernel.
using boundary_t = uint8_t;

/// Index into Lattice::obstacles. types::FLUID marks a fluid node.
/// See boundaries/types.hpp for the sentinel and for solid_mask_t.
using obstacle_id_t = uint16_t;

} // namespace lbm::types
