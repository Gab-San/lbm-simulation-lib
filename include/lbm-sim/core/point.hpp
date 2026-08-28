#ifndef __CORE_POINT_HPP
#define __CORE_POINT_HPP

#include "lbm-sim/cuda/annotations.hpp"
#include "lbm-sim/types/base.hpp"

// C++ STANDARD LIB
#include <iostream>

namespace lbm {

namespace assertion {
template <unsigned short int dim> constexpr bool always_false = false;
}

namespace utils {

template <typename T, types::dim_t dim> struct Point;

template <typename T> struct Point<T, 2> {
  const T x, y;

  LBM_HD_FUNC constexpr Point(const T x_, const T y_) : x(x_), y(y_) {}

  template <typename U>
  LBM_HD_FUNC explicit Point(const Point<U, 2> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)) {}

  template <typename Container>
  constexpr Point(Container c) : Point(c[0], c[1]) {}

  ~Point() = default;
};

template <typename T> struct Point<T, 3> {
  const T x, y, z;

  LBM_HD_FUNC constexpr Point(const T x_, const T y_, const T z_)
      : x(x_), y(y_), z(z_){};

  template <typename U>
  LBM_HD_FUNC explicit Point(const Point<U, 3> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)),
        z(static_cast<T>(other.z)) {}

  template <typename Container>
  constexpr Point(Container c) : Point(c[0], c[1], c[2]) {}

  ~Point() = default;
};

// ---- OPERATOR OVERLOADING ----
template <typename T, types::dim_t dim>
LBM_HD_FUNC inline bool operator==(const Point<T, dim> &lhs,
                                   const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  } else {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }
}

template <typename T, typename K, types::dim_t dim>
LBM_HD_FUNC inline Point<std::common_type_t<T, K>, dim>
operator+(const Point<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  // WARN: Geometrically sum and difference of points does not make sense.
  // All instances where this occurs should be either changed for point by point
  // difference, either should be changed to Point +- Vector.
  if constexpr (dim == 2) {
    return Point<R, 2>(lhs.x + rhs.x, lhs.y + rhs.y);
  } else {
    return Point<R, 3>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }
}

template <typename T, types::dim_t dim>
LBM_HD_FUNC inline Point<T, dim> operator-(const Point<T, dim> &lhs,
                                           const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return Point<T, 2>(lhs.x - rhs.x, lhs.y - rhs.y);
  } else {
    return Point<T, 3>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
}

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
