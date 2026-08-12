#ifndef __CORE_POINT_HPP
#define __CORE_POINT_HPP

// C++ STANDARD LIB
#include <iostream>

namespace CollisionDetection {
namespace utils {

template <typename T, unsigned int dim> struct Point;

template <typename T> struct Point<T, 2> {
  const T x, y;

  Point(const T x_, const T y_) : x(x_), y(y_) {}
  ~Point() = default;

  T measure() const { return x * y; }
};

template <typename T> struct Point<T, 3> {
  const T x, y, z;

  Point(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {};
  ~Point() = default;

  T measure() const { return x * y * z; }
};

// ---- OPERATOR OVERLOADING ----
template <typename T, unsigned int dim>
inline bool operator==(const Point<T, dim> &lhs, const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  } else {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }
}

template <typename T, unsigned int dim>
inline Point<T, dim> operator+(const Point<T, dim> &lhs,
                               const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return Point<T, 2>(lhs.x + rhs.x, lhs.y + rhs.y);
  } else {
    return Point<T, 3>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }
}

template <typename T, unsigned int dim>
inline Point<T, dim> operator-(const Point<T, dim> &lhs,
                               const Point<T, dim> &rhs) {
  if constexpr (dim == 2) {
    return Point<T, 2>(lhs.x - rhs.x, lhs.y - rhs.y);
  } else {
    return Point<T, 3>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
}

template <typename T, unsigned int dim>
inline std::ostream &operator<<(std::ostream &out, const Point<T, dim> &p) {
  if constexpr (dim == 2) {
    return out << "Point(" << p.x << "," << p.y << ")";
  } else {
    return out << "Point(" << p.x << "," << p.y << "," << p.z << ")";
  }
}

} // namespace utils
} // namespace CollisionDetection

#endif // __CORE_POINT_HPP
