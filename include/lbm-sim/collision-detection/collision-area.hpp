#ifndef __COLLISION_AREA_HPP
#define __COLLISION_AREA_HPP

#include "lbm-sim/collision-detection/shape.hpp"

#include "lbm-sim/types/common.hpp"

// C++ STANDARD LIB
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace lbm {
namespace CollisionDetection {

/// Airfoil is 2D-only and CylindricalShell is 3D-only (both static_assert on
/// dim), so neither is an alternative on the other side: std::visit
/// instantiates *every* alternative, and CollisionArea<dim> is visited by
/// Solid::compute_solid_mask<dim>().
template <unsigned short int dim>
using CollisionShapesT = std::conditional_t<
    dim == 2, std::variant<Segment<2>, Circle<2>, Parallelogram<2>, Airfoil<2>>,
    std::variant<Segment<3>, Circle<3>, Parallelogram<3>, CylindricalShell<3>>>;

template <unsigned short int dim> class CollisionArea {

  types::Coordinate<dim> position;

public:
  const std::vector<CollisionShapesT<dim>> collision_shapes;
  
  CollisionArea(const types::Coordinate<dim> position_,
                std::vector<CollisionShapesT<dim>> collision_shapes_)
      : position(position_), collision_shapes(std::move(collision_shapes_)) {}

  ~CollisionArea() = default;

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

  /// Inclusive integer bounding box of the whole area, in grid coordinates.
  /// An area with no shapes gets an inverted box, which
  /// Solid::compute_solid_mask() reads as "nothing to rasterize".
  AABB<dim> aabb() const {
    if (collision_shapes.empty()) {
      types::Coordinate<dim> lo = position;
      types::Coordinate<dim> hi = position;
      for (types::dim_t d = 0; d < dim; ++d) {
        utils::ops::axis(lo, d) = 1;
        utils::ops::axis(hi, d) = 0;
      }
      return {lo, hi};
    }

    bool first = true;
    AABB<dim> box{position, position};
    for (const auto &sh : collision_shapes) {
      const AABB<dim> local =
          std::visit([](const auto &shape) { return shape.aabb(); }, sh);

      // Shapes are defined relative to `position`; contains() subtracts it, so
      // the box has to add it back.
      const AABB<dim> shifted{position + local.min, position + local.max};

      box = first ? shifted : detail::merge(box, shifted);
      first = false;
    }
    return box;
  }
};

} // namespace CollisionDetection
} // namespace lbm

#endif // _SOLID_HPP
