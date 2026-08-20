#ifndef __CORE_VECTOR_HPP
#define __CORE_VECTOR_HPP

#include "lbm-sim/core/point.hpp"

#include "lbm-sim/backend/cuda-annotations.hpp"

// C++ STANDARD LIB
#include <iostream>
#include <type_traits>
#include <utility>

namespace lbm {
namespace utils {

template <typename T, unsigned short int dim> struct Vector;

template <typename T> struct Vector<T, 2> {
  T dx, dy;

  LBM_HD_FUNC constexpr Vector() : dx(0), dy(0) {}
  LBM_HD_FUNC constexpr Vector(const T &x, const T &y) : dx(x), dy(y) {}
  LBM_HD_FUNC constexpr Vector(const Point<T, 2> &A_, const Point<T, 2> &B_)
      : Vector(B_.x - A_.x, B_.y - A_.y) {}

  LBM_HD_FUNC ~Vector() = default;
};

// FIXME: make add correct dimension implementation

// WARN: template may need to be <typename T, int dim, typename K>

// ---- OPERATOR OVERLOADING ----
template <typename T, typename K>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, 2>
operator+(const Vector<T, 2> &lhs, const Vector<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Vector<R, 2>(lhs.dx + rhs.dx, lhs.dy + rhs.dy);
}

template <typename T, typename K>
LBM_HD_FUNC inline Vector<T, 2> &operator+=(Vector<T, 2> &lhs,
                                            const Vector<K, 2> &rhs) {
  lhs.dx += rhs.dx;
  lhs.dy += rhs.dy;
  return lhs;
}

template <typename T>
LBM_HD_FUNC inline Vector<T, 2> operator-(const Vector<T, 2> &lhs,
                                          const Vector<T, 2> &rhs) {
  return Vector<T, 2>(lhs.dx - rhs.dx, lhs.dy - rhs.dy);
}

template <typename T, typename K>
inline Vector<std::common_type_t<T, K>, 2> operator%(const Vector<T, 2> &lhs,
                                                     const Vector<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Vector<R, 2>(lhs.dx % rhs.dx, lhs.dy % rhs.dy);
}

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const Vector<T, 2> &v) {
  return out << "Vector(" << v.dx << "," << v.dy << ")";
}

// ---- OPERATION BY SCALAR ----
template <typename T, typename K>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, 2>
operator*(const Vector<T, 2> &lhs, const K &scalar) {
  using R = std::common_type_t<T, K>;
  return Vector<R, 2>(lhs.dx * scalar, lhs.dy * scalar);
}

template <typename T, typename K>
LBM_HD_FUNC inline Vector<T, 2> &operator/=(Vector<T, 2> &lhs,
                                            const K &scalar) {
  lhs.dx /= scalar;
  lhs.dy /= scalar;
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
