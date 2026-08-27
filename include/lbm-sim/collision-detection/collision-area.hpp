#ifndef __COLLISION_AREA_HPP
#define __COLLISION_AREA_HPP

#include "lbm-sim/collision-detection/shape.hpp"

#include "lbm-sim/types/common.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <utility>
#include <variant>
#include <vector>

namespace lbm {
namespace CollisionDetection {

template <unsigned short int dim>
using CollisionShapesT =
    std::variant<Segment<dim>, Circle<dim>, Parallelogram<dim>, Airfoil<dim>>;

template <unsigned short int dim> class CollisionArea {

  types::Coordinate<dim> position;
  mutable std::vector<types::Coordinate<dim>> cached_perimeter;

public:
  const std::vector<CollisionShapesT<dim>> collision_shapes;
  CollisionArea(const types::Coordinate<dim> position_,
                std::vector<CollisionShapesT<dim>> collision_shapes_)
      : position(position_), collision_shapes(std::move(collision_shapes_)) {}

  ~CollisionArea() = default;

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    for (const auto &sh : collision_shapes) {
      bool hit = std::visit(
          [&](const auto &shape) {
            return shape.isCollidingWith(point - position);
          },
          sh);

      if (hit)
        return true;
    }
    return false;
  };

  bool contains(const types::Coordinate<dim> &point) const {
    for (const auto &sh : collision_shapes) {
      bool hit = std::visit(
          [&](const auto &shape) { return shape.contains(point - position); },
          sh);

      if (hit)
        return true;
    }
    return false;
  }

  const std::vector<types::Coordinate<dim>> &
  getPerimeter(bool force = false) const {
    if (!cached_perimeter.empty() && !force)
      return cached_perimeter;

    cached_perimeter.clear();
    for (const auto &sh : collision_shapes) {
      std::visit(
          [&](const auto &shape) {
            const auto &shape_perimeter = shape.getPerimeter();
            for (const auto &p : shape_perimeter) {
              const auto candidate = position + p;
              if (std::find(cached_perimeter.begin(), cached_perimeter.end(),
                            candidate) != cached_perimeter.end()) {
                continue;
              }

              cached_perimeter.emplace_back(candidate);
            }
          },
          sh);
    }
    return cached_perimeter;
  }
};

} // namespace CollisionDetection
} // namespace lbm

#endif // _SOLID_HPP
