/**
 * @file point.hpp
 * @brief Point<T, dim>: a location in the domain, and the arithmetic on it.
 *
 * Specialised for @c dim == 2 and @c dim == 3 rather than written generically
 * over an array, so the components keep the names @c x, @c y, @c z that the
 * physics is written in, and so a Point stays a trivially copyable aggregate
 * that a CUDA kernel can take by value.
 *
 * The library uses two instantiations, both aliased in types/common.hpp:
 * - @c types::Coordinate<dim> -- @c Point<int, dim>, a lattice node. Signed,
 *   because link resolution routinely produces a negative component before
 *   deciding what the boundary does with it;
 * - @c types::DimPoint<dim> -- @c Point<std::size_t, dim>, a size.
 *
 * @see vector.hpp for the displacement type, and operators.hpp for the mixed
 *      Point/Vector arithmetic.
 */

#ifndef __CORE_POINT_HPP
#define __CORE_POINT_HPP

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/types/fwd.hpp"

// C++ STANDARD LIB
#include <iostream>
#include <type_traits>

namespace lbm {

namespace assertion {
/// Always @c false, but only after @c dim is substituted: the idiom that
/// makes a @c static_assert in the dead branch of an @c if @c constexpr fire
/// on instantiation rather than on parsing.
template <unsigned short int dim> constexpr bool always_false = false;
} // namespace assertion

namespace utils {

/// @brief A location in a 2D domain.
/// @tparam T Component type: @c int for a node, @c std::size_t for a size.
template <typename T> struct Point<T, 2> {
  T x, y; ///< Components, in domain order.

  /// @brief Builds a point from its components.
  LBM_HD_FUNC constexpr Point(const T x_, const T y_) : x(x_), y(y_) {}

  /// @brief Leaves the components uninitialised, so the type stays trivial.
  Point() = default;

  /// @brief Explicit component-wise conversion between component types,
  ///        e.g. DimPoint -> Coordinate. Explicit because it may narrow.
  template <typename U>
  LBM_HD_FUNC explicit Point(const Point<U, 2> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)) {}

  /// @brief Builds from anything indexable, which is how a
  ///        @c std::array<uint64_t,2> read out of a config file becomes a
  ///        point without an intermediate.
  /// @warning Not @c explicit and not constrained, so it is a greedy
  ///          converting constructor: a type with @c operator[] converts
  ///          silently.
  template <typename Container>
  constexpr Point(Container c) : Point(c[0], c[1]) {}

  ~Point() = default;
};

/// @brief A location in a 3D domain.
/// @tparam T Component type: @c int for a node, @c std::size_t for a size.
template <typename T> struct Point<T, 3> {
  T x, y, z; ///< Components, in domain order.

  /// @brief Builds a point from its components.
  LBM_HD_FUNC constexpr Point(const T x_, const T y_, const T z_)
      : x(x_), y(y_), z(z_) {}

  /// @brief Leaves the components uninitialised, so the type stays trivial.
  Point() = default;

  /// @brief Explicit component-wise conversion between component types.
  template <typename U>
  LBM_HD_FUNC explicit Point(const Point<U, 3> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)),
        z(static_cast<T>(other.z)) {}

  /// @brief Builds from anything indexable. Same caveat as the 2D overload.
  template <typename Container>
  constexpr Point(Container c) : Point(c[0], c[1], c[2]) {}

  ~Point() = default;
};

// ---- OPERATOR OVERLOADING ----

/// @brief Component-wise equality.
template <typename T, types::dim_t dim>
LBM_HD_FUNC inline bool operator==(const Point<T, dim> &lhs,
                                   const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  } else {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }
}

/**
 * @brief Component-wise sum of two points.
 *
 * @warning Geometrically meaningless: adding two locations is not a
 *          location. It survives only because a few call sites still use it
 *          as a shorthand for translation; each of those should become
 *          either a point-minus-point difference or a Point +/- Vector.
 *          Do not add new uses.
 */
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator+(const Point<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, 2>(lhs.x + rhs.x, lhs.y + rhs.y);
  } else {
    return Point<R, 3>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }
}

/// @brief Component-wise difference of two points.
/// @note Unlike the sum above this one is meaningful, though it returns a
///       Point where a Vector would be the honest type.
template <typename T, types::dim_t dim>
LBM_HD_FUNC inline Point<T, dim> operator-(const Point<T, dim> &lhs,
                                           const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return Point<T, 2>(lhs.x - rhs.x, lhs.y - rhs.y);
  } else {
    return Point<T, 3>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
}

/// @brief Component-wise modulo, the coordinate wrap behind the periodic
///        boundary conditions.
/// @warning Plain @c %, so a negative component yields a negative result.
///          Callers that may go below zero add the extent first.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator%(const Point<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, dim>(lhs.x % rhs.x, lhs.y % rhs.y);
  } else {
    return Point<R, dim>(lhs.x % rhs.x, lhs.y % rhs.y, lhs.z % rhs.z);
  }
}

/// @brief Stream insertion, for logging: `Point(x,y)` / `Point(x,y,z)`.
///
/// Not @c LBM_HD_FUNC: it is host-only by construction. If a device
/// translation unit ever fails on it, fence it behind the
/// @c __CUDA_ARCH__ guard left commented below.
// NOTE: If compiler returns an error uncomment below
//
// #ifndef __CUDA_ARCH__
template <typename T, types::dim_t dim>
inline std::ostream &operator<<(std::ostream &out, const Point<T, dim> &p) {
  if constexpr (dim == 2) {
    return out << "Point(" << p.x << "," << p.y << ")";
  } else {
    return out << "Point(" << p.x << "," << p.y << "," << p.z << ")";
  }
}
// #endif

} // namespace utils
} // namespace lbm

#endif // __CORE_POINT_HPP
