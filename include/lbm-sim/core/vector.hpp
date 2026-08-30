/**
 * @file vector.hpp
 * @brief Vector<T, dim>: a displacement, and the arithmetic on it.
 *
 * The counterpart of Point: a Point is where something is, a Vector is how
 * far it moves. Components are named @c dx, @c dy, @c dz to keep the two
 * apart at a glance and to make a mix-up a compile error rather than a
 * plausible-looking result.
 *
 * Two instantiations dominate:
 * - @c types::Direction<dim> -- @c Vector<int, dim>, a discrete lattice
 *   velocity @f$ \mathbf{c}_i @f$;
 * - @c Vector<double, dim> -- a macroscopic velocity.
 *
 * Their product is exactly the mixed-type case the @c static_assert at the
 * bottom of this file guards: an integer direction times a @c double field
 * must promote to @c double, not truncate to @c int. Every binary operator
 * therefore returns @c std::common_type_t of its operands.
 *
 * @see point.hpp for locations, and operators.hpp for the mixed Point/Vector
 *      arithmetic and the free algebra (dot, cross).
 */

#ifndef __CORE_VECTOR_HPP
#define __CORE_VECTOR_HPP

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/types/fwd.hpp"

// C++ STANDARD LIB
#include <iostream>
#include <type_traits>
#include <utility>

namespace lbm {
namespace utils {

/// @brief A displacement in a 2D domain.
/// @tparam T Component type: @c int for a lattice direction, @c double for a
///           macroscopic velocity.
template <typename T> struct Vector<T, 2> {
  T dx, dy; ///< Components, in domain order.

  /// @brief Builds a vector from its components.
  LBM_HD_FUNC constexpr Vector(const T &x, const T &y) : dx(x), dy(y) {}

  /// @brief The displacement from @p A_ to @p B_, i.e. @c B_ - A_.
  LBM_HD_FUNC constexpr Vector(const Point<T, 2> &A_, const Point<T, 2> &B_)
      : Vector(B_.x - A_.x, B_.y - A_.y) {}

  /// @brief Builds from anything indexable, which is how the velocity read
  ///        out of a config file becomes a vector without an intermediate.
  /// @warning Not @c explicit and not constrained: a type with
  ///          @c operator[] converts silently.
  template <typename Container>
  constexpr Vector(Container c) : Vector(c[0], c[1]) {}

  /// @brief The zero vector. Unlike Point, a default-constructed Vector *is*
  ///        initialised: it is accumulated into in the moment reduction.
  LBM_HD_FUNC constexpr Vector() : Vector(0, 0) {}

  ~Vector() = default;
};

/// @brief A displacement in a 3D domain.
/// @tparam T Component type: @c int for a lattice direction, @c double for a
///           macroscopic velocity.
template <typename T> struct Vector<T, 3> {
  T dx, dy, dz; ///< Components, in domain order.

  /// @brief Builds a vector from its components.
  LBM_HD_FUNC constexpr Vector(const T &x, const T &y, const T &z)
      : dx(x), dy(y), dz(z) {}

  /// @brief The zero vector.
  LBM_HD_FUNC constexpr Vector() : Vector(0, 0, 0) {}

  /// @brief The displacement from @p A_ to @p B_, i.e. @c B_ - A_.
  LBM_HD_FUNC constexpr Vector(const Point<T, 3> &A_, const Point<T, 3> &B_)
      : Vector(B_.x - A_.x, B_.y - A_.y, B_.z - A_.z) {}

  /// @brief Builds from anything indexable. Same caveat as the 2D overload.
  template <typename Container>
  constexpr Vector(Container c) : Vector(c[0], c[1], c[2]) {}

  ~Vector() = default;
};

// ---- OPERATOR OVERLOADING ----

/// @brief Component-wise sum, promoting to the common component type.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, dim>
operator+(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx + rhs.dx, lhs.dy + rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx + rhs.dx, lhs.dy + rhs.dy, lhs.dz + rhs.dz);
  }
}

/// @brief In-place sum. This is the accumulator of the momentum reduction
///        @f$ \rho\mathbf{u} = \sum_i f_i \mathbf{c}_i @f$.
/// @note The result stays in @c T, so accumulating a @c double into a
///       @c Vector<int> truncates. All current call sites accumulate into a
///       @c Vector<double>.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<T, dim> &operator+=(Vector<T, dim> &lhs,
                                              const Vector<K, dim> &rhs) {
  lhs.dx += rhs.dx;
  lhs.dy += rhs.dy;
  if constexpr (dim == 3) {
    lhs.dz += rhs.dz;
  }
  return lhs;
}

/// @brief Component-wise difference, promoting to the common component type.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, dim>
operator-(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx - rhs.dx, lhs.dy - rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx - rhs.dx, lhs.dy - rhs.dy, lhs.dz - rhs.dz);
  }
}

/**
 * @brief Declared as a modulo, but its body computes @c lhs - rhs.
 *
 * @warning This is a defect, not a convention: the implementation is a
 *          verbatim copy of @c operator-. Nothing in the library calls it
 *          today, which is why it has gone unnoticed. Fix it or delete it
 *          before using it -- do not write code against the current
 *          behaviour.
 */
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, dim>
operator%(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx - rhs.dx, lhs.dy - rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx - rhs.dx, lhs.dy - rhs.dy, lhs.dz - rhs.dz);
  }
}

/// @brief Stream insertion, for logging: `Vector(dx,dy)` / `Vector(dx,dy,dz)`.
template <typename T, types::dim_t dim>
inline std::ostream &operator<<(std::ostream &out, const Vector<T, dim> &v) {
  if constexpr (dim == 2) {
    return out << "Vector(" << v.dx << "," << v.dy << ")";
  } else {
    return out << "Vector(" << v.dx << "," << v.dy << "," << v.dz << ")";
  }
}

// ---- OPERATION BY SCALAR ----

/// @brief Scaling. The @f$ f_i \mathbf{c}_i @f$ of the momentum sum: an
///        @c int direction times a @c double population, promoted.
/// @note Only this operand order exists; @c scalar * vector does not compile.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, dim>
operator*(const Vector<T, dim> &lhs, const K &scalar) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx * scalar, lhs.dy * scalar);
  } else {
    return Vector<R, 3>(lhs.dx * scalar, lhs.dy * scalar, lhs.dz * scalar);
  }
}

/// @brief In-place division by a scalar, the @f$ \mathbf{u} = \rho\mathbf{u}
///        / \rho @f$ that closes the moment computation.
/// @warning No check on @p scalar: a node whose density has collapsed to
///          zero produces infinities here rather than an error.
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Vector<T, dim> &operator/=(Vector<T, dim> &lhs,
                                              const K &scalar) {
  lhs.dx /= scalar;
  lhs.dy /= scalar;
  if constexpr (dim == 3) {
    lhs.dz /= scalar;
  }
  return lhs;
}

// ---- COMPILE-TIME SANITY CHECKS ----
// Mixed-type scalar multiplication must promote, not truncate.
static_assert(std::is_same_v<decltype(std::declval<Vector<int, 2>>() * 1.0),
                             Vector<double, 2>>,
              "Vector<int,2> * double must promote to Vector<double,2>");

static_assert(
    std::is_same_v<decltype(std::declval<Vector<int, 2>>() +
                            std::declval<Vector<double, 2>>()),
                   Vector<double, 2>>,
    "Vector<int,2> + Vector<double,2> must promote to Vector<double,2>");

} // namespace utils
} // namespace lbm

#endif // __CORE_VECTOR_HPP
