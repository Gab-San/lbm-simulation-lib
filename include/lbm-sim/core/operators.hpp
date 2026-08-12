#ifndef __CORE_OPEATORS_HPP
#define __CORE_OPEATORS_HPP

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <type_traits>

namespace CollisionDetection {
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

} // namespace utils
} // namespace CollisionDetection

#endif // __CORE_OPEATORS_HPP
