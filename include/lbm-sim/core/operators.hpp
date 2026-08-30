/**
 * @file operators.hpp
 * @brief Mixed Point/Vector arithmetic and the small vector algebra the
 *        solver needs.
 *
 * Point and Vector are kept apart on purpose: a Point is a location, a
 * Vector is a displacement, and the combinations defined here are the
 * meaningful ones -- point +/- vector gives a point.
 *
 * (point.hpp additionally defines point + point, which is geometrically
 * meaningless and carries a warning saying so; it exists only until its
 * remaining call sites are rewritten.)
 *
 * Every result type goes through @c std::common_type_t, so the recurring
 * mixed case of an integer lattice direction times a @c double field
 * promotes to @c double instead of truncating. A @c static_assert below
 * pins that behaviour down.
 *
 * All functions are @c LBM_HD_FUNC and header-only: they are meant to
 * disappear at @c -O2, on both the host and the device.
 */

#ifndef __CORE_OPEATORS_HPP
#define __CORE_OPEATORS_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/backend/cuda/annotations.hpp"

// C++ STANDARD LIB
#include <type_traits>

namespace lbm::utils {

/// @name Point/Vector arithmetic
/// Translating a point by a displacement, in either operand order, and the
/// two subtractions. The result is always a Point.
/// @{

/// @brief Translates @p lhs by @p rhs.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator+(const Point<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, dim>(lhs.x + rhs.dx, lhs.y + rhs.dy);
  } else {
    return Point<R, dim>(lhs.x + rhs.dx, lhs.y + rhs.dy, lhs.z + rhs.dz);
  }
}

/// @brief Translates @p rhs by @p lhs.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator+(const Vector<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, dim>(lhs.dx + rhs.x, lhs.dy + rhs.y);
  } else {
    return Point<R, dim>(lhs.dx + rhs.x, lhs.dy + rhs.y, lhs.dz + rhs.z);
  }
}

/// @brief Translates @p lhs by @c -rhs.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator-(const Point<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, dim>(lhs.x - rhs.dx, lhs.y - rhs.dy);
  } else {
    return Point<R, dim>(lhs.x - rhs.dx, lhs.y - rhs.dy, lhs.z - rhs.dz);
  }
}

/// @brief Component-wise @c lhs - rhs, kept for symmetry with the overload
///        above.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator-(const Vector<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, dim>(lhs.dx - rhs.x, lhs.dy - rhs.y);
  } else {
    return Point<R, dim>(lhs.dx - rhs.x, lhs.dy - rhs.y, lhs.dz - rhs.z);
  }
}

/// @}

/// Free vector algebra, kept out of the operator overload set so the call
/// sites read as the formulas they implement.
namespace ops {

/**
 * @brief Scalar product.
 *
 * The workhorse of the equilibrium: every @f$ \mathbf{c}_i\cdot\mathbf{u} @f$
 * in the collision kernels goes through here, with an integer @p lhs and a
 * @c double @p rhs.
 */
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline std::common_type_t<T, K> dot(const Vector<T, dim> &lhs,
                                                const Vector<K, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy;
  } else {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy + lhs.dz * rhs.dz;
  }
}

/// @brief 2D cross product: the scalar @c z component of the 3D one.
/// Its sign is the orientation of the pair, which is what the point-in-
/// polygon tests in the shape code use.
template <typename T, typename K>
LBM_HD_FUNC inline std::common_type_t<T, K> cross(const Vector<T, 2> &lhs,
                                                  const Vector<K, 2> &rhs) {
  return lhs.dx * rhs.dy - lhs.dy * rhs.dx;
}

/// @brief 3D cross product.
template <typename T, typename K>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, 3>
cross(const Vector<T, 3> &lhs, const Vector<K, 3> &rhs) {
  using R = std::common_type_t<T, K>;
  return Vector(lhs.dy * rhs.dz - lhs.dz * rhs.dy,
                lhs.dz * rhs.dx - lhs.dx * rhs.dz,
                lhs.dx * rhs.dy - lhs.dy * rhs.dx);
}

/**
 * @brief Product of the components: an area in 2D, a volume in 3D.
 *
 * Used to size a field from a DimPoint before a Grid exists -- for instance
 * in Solid::compute_solid_mask(). Once there is a Grid, Grid::getArea() says
 * the same thing.
 *
 * @tparam dim Must be 2 or 3; anything else is a @c static_assert.
 */
template <types::dim_t dim>
LBM_HD_FUNC inline std::size_t measure(const types::DimPoint<dim> p) {
  if constexpr (dim == 2) {
    return p.x * p.y;
  } else if constexpr (dim == 3) {
    return p.x * p.y * p.z;
  } else {
    static_assert(assertion::always_false<dim>,
                  "DimPoint::measure() with dim > 3!");
  }
}

static_assert(
    std::is_same_v<decltype(dot(std::declval<Vector<int, 2>>(),
                                std::declval<Vector<double, 2>>())),
                   double>,
    "dot(Vector<int,2>, Vector<double,2>) must return double, not int");

/**
 * @brief Component @p a of a point, by index rather than by name.
 *
 * @param p Point to index into.
 * @param a Axis: 0 for x, 1 for y, 2 for z.
 * @return Reference to the component, so it can be assigned through.
 *
 * This is what lets the boundary and bounding-box code loop over the axes
 * instead of writing the 2D and 3D cases out twice.
 *
 * @warning No range check: @p a beyond @c dim-1 returns the last component
 *          rather than failing.
 */
template <typename T, types::dim_t dim>
LBM_HD_FUNC inline T &axis(Point<T, dim> &p, const types::dim_t a) {
  if constexpr (dim == 2) {
    return a == 0 ? p.x : p.y;
  } else {
    return a == 0 ? p.x : (a == 1 ? p.y : p.z);
  }
}

/// @brief Read-only overload of axis().
template <typename T, types::dim_t dim>
LBM_HD_FUNC inline const T &axis(const Point<T, dim> &p, const types::dim_t a) {
  if constexpr (dim == 2) {
    return a == 0 ? p.x : p.y;
  } else {
    return a == 0 ? p.x : (a == 1 ? p.y : p.z);
  }
}

} // namespace ops

} // namespace lbm::utils

#endif // __CORE_OPEATORS_HPP
