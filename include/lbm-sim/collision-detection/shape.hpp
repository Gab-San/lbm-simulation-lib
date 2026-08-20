#ifndef __SHAPE_HPP
#define __SHAPE_HPP

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/collision-detection/algorithms/collision.hpp"
#include "lbm-sim/collision-detection/algorithms/rasterization.hpp"

// C++ STANDARD LIB
#include <stdexcept>
#include <vector>

// FIXME: Should be extendible to 3D

namespace lbm {
namespace CollisionDetection {

template <unsigned short int dim, class DerivedShape> class Shape {
public:
  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    return static_cast<DerivedShape *>(this)->isCollidingWith(point);
  };

  bool contains(const types::Coordinate<dim> &point) const {
    return static_cast<DerivedShape *>(this)->contains(point);
  };

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    return static_cast<DerivedShape *>(this)->getPerimeter();
  };

protected:
  Shape() = default;
};

template <unsigned short int dim>
class Segment : public Shape<dim, Segment<dim>> {
  const types::Coordinate<dim> A, B;
  mutable std::vector<types::Coordinate<dim>> cached_perimeter;
  const double dist_x, dist_y;

public:
  Segment(const types::Coordinate<dim> A_, const types::Coordinate<dim> B_)
      : A(A_), B(B_), dist_x(static_cast<double>(B.x - A.x)),
        dist_y(static_cast<double>(B.y - A.y)) {}

  ~Segment() = default;

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    // const auto &perimeter = getPerimeter();
    // return std::find(perimeter.begin(), perimeter.end(), point) !=
    // perimeter.end();

    return algorithms::bounding_box_check(A, B, point) &&
           algorithms::brasenham_collision(A, B, point);
  }

  bool contains(const types::Coordinate<dim> &point) const {
    return isCollidingWith(point);
  }

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    if (!cached_perimeter.empty()) {
      return cached_perimeter;
    }

    const utils::Vector<int, dim> AB(A, B);

    // FIXME: This needs to be extended to 3D case
    if (AB.dx < 0) {
      cached_perimeter = algorithms::brasenham_rasterisation<2>(B, A);
    } else {
      cached_perimeter = algorithms::brasenham_rasterisation<2>(A, B);
    }

    return cached_perimeter;
  }
};

template <unsigned short int dim>
class Circle : public Shape<dim, Circle<dim>> {
  const types::Coordinate<dim> center;
  const unsigned int radius;
  mutable std::vector<types::Coordinate<dim>> cached_perimeter;

public:
  Circle(const types::Coordinate<dim> center_, const unsigned int radius_)
      : center(center_), radius(radius_) {}

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    if constexpr (dim == 2) {
      return checkCollisionWith(point);
    } else {
      return checkCollisionWith(point);
    }
  }

  bool contains(const types::Coordinate<dim> &point) const {
    const int center_point_dist_x = center.x - point.x;
    const unsigned int dist_x = center_point_dist_x * center_point_dist_x;
    const int center_point_dist_y = center.y - point.y;
    const unsigned int dist_y = center_point_dist_y * center_point_dist_y;

    return dist_x + dist_y <= radius * radius;
  }

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    if (!cached_perimeter.empty()) {
      return cached_perimeter;
    }

    cached_perimeter.clear();

    int x = -radius;
    int y = 0;
    cached_perimeter.push_back(center + types::Coordinate<dim>(x, y));
    cached_perimeter.push_back(center + types::Coordinate<dim>(-x, y));

    cached_perimeter.push_back(center + types::Coordinate<dim>(0, radius));
    cached_perimeter.push_back(center + types::Coordinate<dim>(0, -radius));

    while (x < static_cast<int>(radius)) {
      if (x < 0) {
        if (x * x + (y + 1) * (y + 1) > radius * radius) {
          x++;
        } else {
          y++;
        }
      } else {
        if ((x + 1) * (x + 1) + y * y > radius * radius) {
          y--;
        } else {
          x++;
        }
      }

      if (y == -y) {
        continue;
      }

      cached_perimeter.push_back(center + types::Coordinate<dim>(x, y));
      cached_perimeter.push_back(center + types::Coordinate<dim>(x, -y));
    }

    return cached_perimeter;
  }

private:
  bool checkCollisionWith(const types::Coordinate<2> &point) const {
    const int center_point_dist_x = center.x - point.x;
    const unsigned int dist_x = center_point_dist_x * center_point_dist_x;
    const int center_point_dist_y = center.y - point.y;
    const unsigned int dist_y = center_point_dist_y * center_point_dist_y;

    return dist_x + dist_y == radius * radius;
  }

  bool checkCollisionWith(const types::Coordinate<3> &point) const {
    throw std::runtime_error("Circle::checkCollisionWith(utils::Point<3>&) : "
                             "Collision Detection in 3D not yet implemented!");
  }
};

template <unsigned short int dim>
class Parallelogram : public Shape<dim, Parallelogram<dim>> {
  const types::Coordinate<dim> A, B, C, D;
  mutable std::vector<types::Coordinate<dim>> cached_perimeter;

public:
  Parallelogram(const types::Coordinate<dim> A_,
                const types::Coordinate<dim> B_,
                const types::Coordinate<dim> C_,
                const types::Coordinate<dim> D_)
      : A(A_), B(B_), C(C_), D(D_) {
    utils::Vector<int, dim> AB(A_, B_), BC(B_, C_), CD(C_, D_), DA(D_, A_);
    if (cross(AB, CD) != 0 || cross(BC, DA) != 0) {
      throw std::invalid_argument(
          "Parallelogram : opposite sides not parallel");
    }
  }

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    Segment<dim> sideAB(A, B);
    Segment<dim> sideBC(B, C);
    Segment<dim> sideCD(C, D);
    Segment<dim> sideDA(D, A);

    bool hit = sideAB.isCollidingWith(point) || sideBC.isCollidingWith(point) ||
               sideCD.isCollidingWith(point) || sideDA.isCollidingWith(point);

    return hit;
  }

  bool contains(const types::Coordinate<dim> &point) const {
    const utils::Vector<int, dim> AB(A, B), AD(A, D), AP(A, point);

    const int dotAB_AB = dot(AB, AB);
    const int dotAP_AB = dot(AP, AB);
    const int dotAD_AD = dot(AD, AD);
    const int dotAP_AD = dot(AP, AD);

    // The dot product can be considered the module of the projection
    // of vector AP onto the target vector (AB, AD).
    // The projection of a vector onto itself (if normalized) is equal to 1.
    // The projection of any other vector will be between [-1, 1).
    //
    // In this case we check that the projection of the vector is positive and
    // it is less then the projection of the target vector onto itself (if
    // normalized this would correspond to checking that the projection lies in
    // [0,1]).
    return dotAP_AB >= 0 && dotAP_AB <= dotAB_AB && dotAP_AD >= 0 &&
           dotAP_AD <= dotAD_AD;
  }

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    if (!cached_perimeter.empty()) {
      return cached_perimeter;
    }

    cached_perimeter.clear();
    Segment<dim> sideAB(A, B);
    Segment<dim> sideBC(B, C);
    Segment<dim> sideCD(C, D);
    Segment<dim> sideDA(D, A);

    for (const auto &p : sideAB.getPerimeter()) {
      cached_perimeter.push_back(p);
    }
    for (const auto &p : sideBC.getPerimeter()) {
      if (p == B)
        continue;
      cached_perimeter.push_back(p);
    }
    for (const auto &p : sideCD.getPerimeter()) {
      if (p == C)
        continue;
      cached_perimeter.push_back(p);
    }
    for (const auto &p : sideDA.getPerimeter()) {
      if (p == D || p == A)
        continue;
      cached_perimeter.push_back(p);
    }

    return cached_perimeter;
  }
};

} // namespace CollisionDetection
} // namespace lbm

#endif // __SHAPE_HPP
