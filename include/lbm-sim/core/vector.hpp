#ifndef __CORE_VECTOR_HPP
#define __CORE_VECTOR_HPP

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/types/base.hpp"

// C++ STANDARD LIB
#include <iostream>
#include <type_traits>
#include <utility>

namespace lbm {
namespace utils {

template <typename T, types::dim_t dim> struct Vector;

template <typename T> struct Vector<T, 2> {
  T dx, dy;

  constexpr Vector(const T &x, const T &y) : dx(x), dy(y) {}

  constexpr Vector() : Vector(0, 0) {}
  constexpr Vector(const Point<T, 2> &A_, const Point<T, 2> &B_)
      : Vector(B_.x - A_.x, B_.y - A_.y) {}

  ~Vector() = default;
};

template <typename T> struct Vector<T, 3> {
  T dx, dy, dz;

  constexpr Vector(const T &x, const T &y, const T &z) : dx(x), dy(y), dz(z) {}

  constexpr Vector() : Vector(0, 0, 0) {}
  constexpr Vector(const Point<T, 3> &A_, const Point<T, 3> &B_)
      : Vector(B_.x - A_.x, B_.y - A_.y, B_.z - A_.z) {}

  ~Vector() = default;
};

// ---- OPERATOR OVERLOADING ----
template <typename T, typename K, types::dim_t dim>
inline Vector<std::common_type_t<T, K>, dim>
operator+(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx + rhs.dx, lhs.dy + rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx + rhs.dx, lhs.dy + rhs.dy, lhs.dz + rhs.dz);
  }
}

template <typename T, typename K, types::dim_t dim>
inline Vector<T, dim> &operator+=(Vector<T, dim> &lhs,
                                  const Vector<K, dim> &rhs) {
  lhs.dx += rhs.dx;
  lhs.dy += rhs.dy;
  if constexpr (dim == 3) {
    lhs.dz += rhs.dz;
  }
  return lhs;
}

template <typename T, typename K, types::dim_t dim>
inline Vector<std::common_type<T, K>, dim>
operator-(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx - rhs.dx, lhs.dy - rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx - rhs.dx, lhs.dy - rhs.dy, lhs.dz - rhs.dz);
  }
}

template <typename T, typename K, types::dim_t dim>
inline Vector<std::common_type_t<T, K>, dim>
operator%(const Vector<T, dim> &lhs, const Vector<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx % rhs.dx, lhs.dy % rhs.dy);
  } else {
    return Vector<R, 3>(lhs.dx % rhs.dx, lhs.dy % rhs.dy, lhs.dz % rhs.dz);
  }
}

template <typename T, types::dim_t dim>
inline std::ostream &operator<<(std::ostream &out, const Vector<T, dim> &v) {
  if constexpr (dim == 2) {
    return out << "Vector(" << v.dx << "," << v.dy << ")";
  } else {
    return out << "Vector(" << v.dx << "," << v.dy << v.dz << ")";
  }
}

// ---- OPERATION BY SCALAR ----
template <typename T, typename K, types::dim_t dim>
inline Vector<std::common_type_t<T, K>, dim>
operator*(const Vector<T, dim> &lhs, const K &scalar) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Vector<R, 2>(lhs.dx * scalar, lhs.dy * scalar);
  } else {
    return Vector<R, 3>(lhs.dx * scalar, lhs.dy * scalar, lhs.dz * scalar);
  }
}

template <typename T, typename K, types::dim_t dim>
inline Vector<T, dim> &operator/=(Vector<T, dim> &lhs, const K &scalar) {
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
