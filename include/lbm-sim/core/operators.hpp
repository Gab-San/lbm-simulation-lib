#ifndef __CORE_OPEATORS_HPP
#define __CORE_OPEATORS_HPP

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <type_traits>

namespace lbm {
namespace utils {

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 2> operator+(const Point<T, 2> &lhs,
                                                    const Vector<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 2>(lhs.x + rhs.dx, lhs.y + rhs.dy);
}

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 2> operator+(const Vector<T, 2> &lhs,
                                                    const Point<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 2>(lhs.x + rhs.dx, lhs.y + rhs.dy);
}

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 2> operator-(const Point<T, 2> &lhs,
                                                    const Vector<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 2>(lhs.x - rhs.dx, lhs.y - rhs.dy);
}

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 2> operator-(const Vector<T, 2> &lhs,
                                                    const Point<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 2>(lhs.x - rhs.dx, lhs.y - rhs.dy);
}

namespace ops {
template <typename T, typename K, unsigned short int dim>
inline std::common_type_t<T, K> dot(const Vector<T, dim> &lhs,
                                    const Vector<K, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy;
  } else {
    static_assert(assertion::always_false<dim>,
                  "dot() : operator not yet implemented for 3D");
  }
}

template <typename T, typename K, unsigned short int dim>
inline std::common_type_t<T, K> cross(const Vector<T, dim> &lhs,
                                      const Vector<K, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.dx * rhs.dy - lhs.dy * rhs.dx;
  } else {
    static_assert(assertion::always_false<dim>,
                  "cross() : operator not yet implemented for 3D");
  }
}

template <unsigned short int dim>
inline std::size_t measure(const types::DimPoint<dim> p) {
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

} // namespace ops

} // namespace utils
} // namespace lbm

#endif // __CORE_OPEATORS_HPP
