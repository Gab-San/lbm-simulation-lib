#ifndef __COLLISION_DETECTION_ALGORITHMS_RASTERIZATION_HPP
#define __COLLISION_DETECTION_ALGORITHMS_RASTERIZATION_HPP

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <vector>

namespace lbm {
namespace CollisionDetection {
namespace algorithms {

// FIXME: Extend to 3D
// FIXME: Extend to all directions
template <unsigned short int dim>
std::vector<types::Coordinate<2>>
brasenham_rasterisation(types::Coordinate<2> A, types::Coordinate<2> B) {
  std::vector<types::Coordinate<2>> rasterized_segment;
  int error = 0;
  const utils::Vector<int, 2> AB(A, B);

  if (AB.dx >= 0 && AB.dy >= 0) { // FIRST QUADRANT

    if (AB.dx >= AB.dy) { // FIRST OCTANT
      // Evaluating 0<=m<=1 : AB.dx > AB.dy
      int y = A.y;

      // if 2 (error + dy) < dx then error += dy;
      // else y++; error += dy - dx;
      //
      // condition: error + dy/dx < 0.5

      for (int x = A.x; x <= B.x; x++) {
        rasterized_segment.push_back({x, y});

        error += AB.dy;
        if ((error << 1) >= AB.dx) {
          y++;
          error -= AB.dx;
        }
      }
    } else { // SECOND OCTANT
      int x = A.x;

      for (int y = A.y; y <= B.y; y++) {
        rasterized_segment.push_back({x, y});

        error += AB.dx;
        if ((error << 1) >= AB.dy) {
          x++;
          error -= AB.dy;
        }
      }
    }
  } else if (AB.dx >= 0 && AB.dy < 0) { // FOURTH QUADRANT
    if (AB.dx >= -AB.dy) {              // EIGHTH OCTANT
      int y = A.y;
      for (int x = A.x; x <= B.x; x++) {
        rasterized_segment.push_back({x, y});

        error += AB.dy;
        if ((error << 1) <= -0.5) {
          y--;
          error += AB.dx;
        }
      }
    } else {
      int x = A.x;

      for (int y = A.y; y >= B.y; y--) {
        rasterized_segment.push_back({x, y});

        error += AB.dy;
        if ((error << 1) >= AB.dx) {
          x++;
          error += AB.dx;
        }
      }
    }
  }

  return rasterized_segment;
}

} // namespace algorithms
} // namespace CollisionDetection
} // namespace lbm

#endif // __COLLISION_DETECTION_ALGORITHMS_RASTERIZATION_HPP
