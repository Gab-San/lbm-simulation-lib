#ifndef __CORE_POINT_HPP
#define __CORE_POINT_HPP

// C++ STANDARD LIB
#include <iostream>

namespace lbm {

namespace assertion {
template <unsigned short int dim> constexpr bool always_false = false;
}

namespace utils {

template <typename T, unsigned short int dim> struct Point;

template <typename T> struct Point<T, 2> {
  const T x, y;

  Point(const T x_, const T y_) : x(x_), y(y_) {}
  template <typename U>
  explicit Point(const Point<U, 2> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)) {}
  ~Point() = default;
};

template <typename T> struct Point<T, 3> {
  const T x, y, z;

  Point(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {};
  template <typename U>
  explicit Point(const Point<U, 3> &other)
      : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)),
        z(static_cast<T>(other.z)) {}
  ~Point() = default;
};

// ---- OPERATOR OVERLOADING ----
template <typename T, unsigned short int dim>
inline bool operator==(const Point<T, dim> &lhs, const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  } else {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }
}

template <typename T, typename K, unsigned short int dim>
inline Point<std::common_type_t<T, K>, dim>
operator+(const Point<T, dim> &lhs, const Point<K, dim> &rhs) {
  using R = std::common_type_t<T, K>;
  if constexpr (dim == 2) {
    return Point<R, 2>(lhs.x + rhs.x, lhs.y + rhs.y);
  } else {
    return Point<R, 3>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }
}

template <typename T, unsigned short int dim>
inline Point<T, dim> operator-(const Point<T, dim> &lhs,
                               const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return Point<T, 2>(lhs.x - rhs.x, lhs.y - rhs.y);
  } else {
    return Point<T, 3>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
}

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 2> operator%(const Point<T, 2> &lhs,
                                                    const Point<K, 2> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 2>(lhs.x % rhs.x, lhs.y % rhs.y);
}

template <typename T, typename K>
inline Point<std::common_type_t<T, K>, 3> operator%(const Point<T, 3> &lhs,
                                                    const Point<K, 3> &rhs) {
  using R = std::common_type_t<T, K>;
  return Point<R, 3>(lhs.x % rhs.x, lhs.y % rhs.y, lhs.z % rhs.z);
}

template <typename T, unsigned short int dim>
inline std::ostream &operator<<(std::ostream &out, const Point<T, dim> &p) {
  if constexpr (dim == 2) {
    return out << "Point(" << p.x << "," << p.y << ")";
  } else {
    return out << "Point(" << p.x << "," << p.y << "," << p.z << ")";
  }
}

} // namespace utils
} // namespace lbm

#endif // __CORE_POINT_HPP
