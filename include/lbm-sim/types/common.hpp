/**
 * @file common.hpp
 * @brief The named Point/Vector instantiations the library actually uses.
 *
 * Four aliases, and the distinctions between them are load-bearing:
 * a size is unsigned, a node coordinate is signed (link resolution goes
 * negative before the boundary decides what to do about it), and a direction
 * is a Vector rather than a Point because it is a displacement.
 *
 * This is the header to include to get the complete types; types/fwd.hpp is
 * the declaration-only counterpart.
 */

#pragma once

// The aliases below name Point/Vector as complete types, so pull in the
// definitions, not just types/fwd.hpp: every user of Coordinate<dim> et al.
// would otherwise have to remember to include core/point.hpp itself. Neither
// header includes this one, so there is no cycle.
#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/types/base.hpp"
#include "lbm-sim/types/fwd.hpp"

#include <cstddef>

namespace lbm::types {

/// A location with continuous coordinates, for geometry that is not snapped
/// to the lattice.
template <dim_t dim> using ValuePoint = utils::Point<double, dim>;

/// Domain extents: a cell count per axis. Unsigned, since it sizes fields.
template <dim_t dim> using DimPoint = utils::Point<std::size_t, dim>;

/// A lattice node. Signed on purpose: adding a direction to a node on the
/// border yields a negative component, which is exactly the case the
/// boundary conditions have to see before it is wrapped or reflected.
template <dim_t dim> using Coordinate = utils::Point<int, dim>;

/// A discrete lattice velocity @f$ \mathbf{c}_i @f$. Components are in
/// @c {-1, 0, 1} for every velocity set in the library.
template <dim_t dim> using Direction = utils::Vector<int, dim>;

} // namespace lbm::types
