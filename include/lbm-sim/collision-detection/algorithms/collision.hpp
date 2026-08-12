#ifndef __COLLISION_DETECTION_ALGORITHMS_COLLISION_HPP
#define __COLLISION_DETECTION_ALGORITHMS_COLLISION_HPP

#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STD LIB
#include <assert.h>

namespace lbm {
namespace CollisionDetection {
namespace algorithms {

bool bounding_box_check(types::Coordinate<2> A, types::Coordinate<2> B,
                        types::Coordinate<2> P) {
  const int maxx = A.x < B.x ? B.x : A.x;
  const int minx = A.x < B.x ? A.x : B.x;
  const int maxy = A.y < B.y ? B.y : A.y;
  const int miny = A.y < B.y ? A.y : B.y;

  return minx <= P.x && P.x <= maxx && miny <= P.y && P.y <= maxy;
}

bool brasenham_collision(types::Coordinate<2> A, types::Coordinate<2> B,
                         types::Coordinate<2> P) {

  const utils::Vector<int, 2> AB(A, B);
  const utils::Vector<int, 2> AP(A, P);

  const int cross_product = utils::ops::cross(AP, AB);
  return 2 * std::abs(cross_product) <
         std::max(std::abs(AB.dx), std::abs(AB.dy));
}
} // namespace algorithms
} // namespace CollisionDetection
} // namespace lbm

#endif // __COLLISION_DETECTION_ALGORITHMS_COLLISION_HPP
