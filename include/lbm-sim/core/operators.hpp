#ifndef __CORE_OPEATORS_HPP
#define __CORE_OPEATORS_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/cuda/annotations.hpp"

// C++ STANDARD LIB
#include <type_traits>

namespace lbm::utils {

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

namespace ops {
template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline std::common_type_t<T, K> dot(const Vector<T, dim> &lhs,
                                                const Vector<K, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy;
  } else {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy + lhs.dz * rhs.dz;
  }
}

template <typename T, typename K>
LBM_HD_FUNC inline std::common_type_t<T, K> cross(const Vector<T, 2> &lhs,
                                                  const Vector<K, 2> &rhs) {
  return lhs.dx * rhs.dy - lhs.dy * rhs.dx;
}

template <typename T, typename K>
LBM_HD_FUNC inline Vector<std::common_type_t<T, K>, 3>
cross(const Vector<T, 3> &lhs, const Vector<K, 3> &rhs) {
  using R = std::common_type_t<T, K>;
  return Vector(lhs.dy * rhs.dz - lhs.dz * rhs.dy,
                lhs.dz * rhs.dx - lhs.dx * rhs.dz,
                lhs.dx * rhs.dy - lhs.dy * rhs.dx);
}

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

template <typename T, types::dim_t dim>
LBM_HD_FUNC inline T &axis(Point<T, dim> &p, const types::dim_t a) {
  if constexpr (dim == 2) {
    return a == 0 ? p.x : p.y;
  } else {
    return a == 0 ? p.x : (a == 1 ? p.y : p.z);
  }
}

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
