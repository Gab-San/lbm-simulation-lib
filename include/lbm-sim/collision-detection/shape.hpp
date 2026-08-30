/**
 * @file shape.hpp
 * @brief The analytic solid shapes, and the CRTP base they share.
 *
 * A shape answers two questions and no others: is this lattice node inside
 * the solid (@c contains()), and which box is worth scanning to find out
 * (@c aabb()). compute_solid_mask() needs nothing else to rasterise it.
 *
 * Dispatch is CRTP rather than virtual, so a shape can be stored by value in
 * the @c std::variant that CollisionArea holds and visited without an
 * indirect call. The price is that the set of shapes is closed at compile
 * time: adding one means adding it to CollisionShapesT<dim> as well.
 *
 * @note @c aabb() is inclusive on both ends and unclamped;
 *       compute_solid_mask() clamps it against the grid. Returning a box
 *       that is too large only costs setup time, one that is too small
 *       silently truncates the body.
 *
 * @see collision-area.hpp for how shapes are grouped and positioned.
 */

#ifndef __SHAPE_HPP
#define __SHAPE_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/collision-detection/algorithms/collision.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lbm {
namespace CollisionDetection {

/// Integer bounding box of a shape. INCLUSIVE on both ends and unclamped:
/// Solid::compute_solid_mask() clamps it against the grid.
template <types::dim_t dim> struct AABB {
  types::Coordinate<dim> min, max;
};

namespace detail {

/// Component-wise hull of two inclusive boxes.
template <types::dim_t dim>
inline AABB<dim> merge(const AABB<dim> &a, const AABB<dim> &b) {
  using utils::ops::axis;

  AABB<dim> out = a;
  for (types::dim_t d = 0; d < dim; ++d) {
    int &lo = axis(out.min, d);
    int &hi = axis(out.max, d);
    if (axis(b.min, d) < lo)
      lo = axis(b.min, d);
    if (axis(b.max, d) > hi)
      hi = axis(b.max, d);
  }
  return out;
}

/// Inclusive box spanning a handful of points.
template <types::dim_t dim>
inline AABB<dim> hull(std::initializer_list<types::Coordinate<dim>> pts) {
  using utils::ops::axis;

  AABB<dim> box{*pts.begin(), *pts.begin()};
  for (const types::Coordinate<dim> &v : pts) {
    for (types::dim_t d = 0; d < dim; ++d) {
      int &lo = axis(box.min, d);
      int &hi = axis(box.max, d);
      if (axis(v, d) < lo)
        lo = axis(v, d);
      if (axis(v, d) > hi)
        hi = axis(v, d);
    }
  }
  return box;
}

} // namespace detail

template <types::dim_t dim, class DerivedShape> class Shape {
public:
  /// True when the lattice node `point` is inside the solid.
  bool contains(const types::Coordinate<dim> &point) const {
    return static_cast<const DerivedShape *>(this)->contains(point);
  };

  /// Inclusive integer bounding box, in the same frame as contains().
  AABB<dim> aabb() const {
    return static_cast<const DerivedShape *>(this)->aabb();
  };

protected:
  Shape() = default;
};

/**
 * Zero-measure shape, and the one documented exception to the
 * "contains() is an exact inequality" rule: nothing is strictly inside a 1-D
 * set, so a rasterized line test stands in for it. This is not licence to
 * reintroduce Bresenham for shapes that do have an interior.
 */
template <types::dim_t dim> class Segment : public Shape<dim, Segment<dim>> {
  const types::Coordinate<dim> A, B;

public:
  Segment(const types::Coordinate<dim> A_, const types::Coordinate<dim> B_)
      : A(A_), B(B_) {}

  ~Segment() = default;

  bool contains(const types::Coordinate<dim> &point) const {
    if constexpr (dim == 2) {
      return algorithms::bounding_box_check(A, B, point) &&
             algorithms::brasenham_collision(A, B, point);
    } else {
      // Reachable only if a Segment<3> ends up in an obstacle list: the visit
      // in CollisionArea::contains() instantiates every alternative, so this
      // has to compile, but the line rasterizer is 2D-only.
      throw std::runtime_error("Segment<3>::contains() : segment collision "
                               "detection in 3D not yet implemented!");
    }
  }

  AABB<dim> aabb() const { return detail::hull<dim>({A, B}); }
};

template <types::dim_t dim> class Circle : public Shape<dim, Circle<dim>> {
  const types::Coordinate<dim> center;
  const unsigned int radius;

public:
  Circle(const types::Coordinate<dim> center_, const unsigned int radius_)
      : center(center_), radius(radius_) {}

  /// Exact integer test: dist^2 <= r^2. No sqrt, and no `==` against r^2 --
  /// essentially no lattice point sits exactly on the circle.
  bool contains(const types::Coordinate<dim> &point) const {
    const types::Direction<dim> d(center, point);
    const int r = static_cast<int>(radius);
    return utils::ops::dot(d, d) <= r * r;
  }

  AABB<dim> aabb() const {
    const int r = static_cast<int>(radius);
    if constexpr (dim == 2) {
      return {{center.x - r, center.y - r}, {center.x + r, center.y + r}};
    } else {
      return {{center.x - r, center.y - r, center.z - r},
              {center.x + r, center.y + r, center.z + r}};
    }
  }
};

template <types::dim_t dim>
class Parallelogram : public Shape<dim, Parallelogram<dim>> {
  const types::Coordinate<dim> A, B, C, D;

public:
  Parallelogram(const types::Coordinate<dim> A_,
                const types::Coordinate<dim> B_,
                const types::Coordinate<dim> C_,
                const types::Coordinate<dim> D_)
      : A(A_), B(B_), C(C_), D(D_) {
    utils::Vector<int, dim> AB(A_, B_), BC(B_, C_), CD(C_, D_), DA(D_, A_);
    if (!is_parallel(AB, CD) || !is_parallel(BC, DA)) {
      throw std::invalid_argument(
          "Parallelogram : opposite sides not parallel");
    }
  }

  bool contains(const types::Coordinate<dim> &point) const {
    const utils::Vector<int, dim> AB(A, B), AD(A, D), AP(A, point);

    const int dotAB_AB = utils::ops::dot(AB, AB);
    const int dotAP_AB = utils::ops::dot(AP, AB);
    const int dotAD_AD = utils::ops::dot(AD, AD);
    const int dotAP_AD = utils::ops::dot(AP, AD);

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

  AABB<dim> aabb() const { return detail::hull<dim>({A, B, C, D}); }

private:
  /// cross() is a scalar in 2D and a vector in 3D, so the zero test has to be
  /// written once per case.
  static bool is_parallel(const utils::Vector<int, dim> &lhs,
                          const utils::Vector<int, dim> &rhs) {
    if constexpr (dim == 2) {
      return utils::ops::cross(lhs, rhs) == 0;
    } else {
      const auto c = utils::ops::cross(lhs, rhs);
      return utils::ops::dot(c, c) == 0;
    }
  }
};

/**
 * Wall of a pipe: the cylindrical shell around the grid-aligned axis A--B.
 * A node is solid when it sits between the two endpoints along that axis and
 * its distance from the axis falls in (inner_radius, outer_radius].
 *
 * This is the COMPLEMENT of a cylinder, and deliberately so: in a duct the
 * fluid is the inside and the solid is everything around it, while
 * Solid::compute_solid_mask() paints exactly what contains() calls solid.
 * Give outer_radius a value large enough to leave the domain -- the AABB is
 * clamped against the grid anyway -- and the mask becomes "every node past
 * inner_radius", i.e. a pipe drilled through the box.
 *
 * Exact integer test, like Circle: dist^2 against the two radii, no sqrt. The
 * accumulator is 64-bit because outer_radius is routinely set well past the
 * grid, and r^2 for a few thousand cells already leaves int territory.
 */
template <types::dim_t dim>
class CylindricalShell : public Shape<dim, CylindricalShell<dim>> {
  static_assert(dim == 3, "CylindricalShell currently supports only 3D");

  const types::Coordinate<dim> A, B;
  const unsigned int inner_radius, outer_radius;
  const types::dim_t axis_id;

  /// The one axis along which the two endpoints differ: that is the pipe
  /// direction. Anything else (a skew axis, or A == B) is rejected here
  /// rather than silently rasterized as something unintended.
  static types::dim_t deduce_axis(const types::Coordinate<dim> &A_,
                                  const types::Coordinate<dim> &B_) {
    types::dim_t found = 0;
    int ndiff = 0;
    for (types::dim_t d = 0; d < dim; ++d) {
      if (utils::ops::axis(A_, d) != utils::ops::axis(B_, d)) {
        found = d;
        ++ndiff;
      }
    }
    if (ndiff != 1) {
      throw std::invalid_argument("CylindricalShell : the two endpoints must "
                                  "differ along exactly one grid axis");
    }
    return found;
  }

public:
  CylindricalShell(const types::Coordinate<dim> A_,
                   const types::Coordinate<dim> B_,
                   const unsigned int inner_radius_,
                   const unsigned int outer_radius_)
      : A(A_), B(B_), inner_radius(inner_radius_), outer_radius(outer_radius_),
        axis_id(deduce_axis(A_, B_)) {
    if (inner_radius_ >= outer_radius_) {
      throw std::invalid_argument(
          "CylindricalShell : inner radius must be smaller than the outer one");
    }
  }

  bool contains(const types::Coordinate<dim> &point) const {
    using utils::ops::axis;

    const int s = axis(point, axis_id);
    const int lo = std::min(axis(A, axis_id), axis(B, axis_id));
    const int hi = std::max(axis(A, axis_id), axis(B, axis_id));
    if (s < lo || s > hi)
      return false;

    // Squared distance from the axis: the endpoints share every coordinate
    // but axis_id, so A carries the centre of the cross section.
    long long r2 = 0;
    for (types::dim_t d = 0; d < dim; ++d) {
      if (d == axis_id)
        continue;
      const long long e = axis(point, d) - axis(A, d);
      r2 += e * e;
    }

    const long long rin = inner_radius;
    const long long rout = outer_radius;
    return r2 > rin * rin && r2 <= rout * rout;
  }

  AABB<dim> aabb() const {
    using utils::ops::axis;

    const int r = static_cast<int>(outer_radius);
    types::Coordinate<dim> lo = A, hi = A;
    for (types::dim_t d = 0; d < dim; ++d) {
      if (d == axis_id) {
        axis(lo, d) = std::min(axis(A, d), axis(B, d));
        axis(hi, d) = std::max(axis(A, d), axis(B, d));
      } else {
        axis(lo, d) = axis(A, d) - r;
        axis(hi, d) = axis(A, d) + r;
      }
    }
    return {lo, hi};
  }
};

template <unsigned short int dim>
class Airfoil : public Shape<dim, Airfoil<dim>> {
  static_assert(dim == 2, "Airfoil currently supports only 2D");

  const types::Coordinate<dim> position; // reference point (leading edge)
  const double chord;                    // chord, in cells
  const double thickness;                // XX/100 (e.g. 0.12 for NACA 0012)
  const double max_camber;               // M/100 (e.g. 0.02 for NACA 2412)
  const double camber_pos;               // P/10 (e.g. 0.4 for NACA 2412)
  const double aoa_rad;                  // angle of attack, in radians

  mutable std::vector<std::pair<double, double>>
      cached_polygon; // normalised (x,y) points of the closed contour

  // Builds the closed contour (upper side, then lower side reversed) in
  // chord-normalised [0,1] coordinates.
  void buildPolygon() const {
    constexpr int N = 100; // resolution: higher = more accurate
    std::vector<std::pair<double, double>> upper, lower;

    for (int i = 0; i <= N; ++i) {
      // Cosine spacing, to cluster points near the leading/trailing edge.
      // M_PI is spelled out for safety, even though cmath should define it.
      double beta = 3.14159265358979323846 * i / N;
      double x = 0.5 * (1 - std::cos(beta));

      double yt = 5.0 * thickness *
                  (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                   0.2843 * x * x * x - 0.1015 * x * x * x * x);

      double yc = 0.0, dyc_dx = 0.0;
      if (max_camber > 0.0 && camber_pos > 0.0) {
        if (x < camber_pos) {
          yc = (max_camber / (camber_pos * camber_pos)) *
               (2 * camber_pos * x - x * x);
          dyc_dx =
              (2 * max_camber / (camber_pos * camber_pos)) * (camber_pos - x);
        } else {
          yc = (max_camber / ((1 - camber_pos) * (1 - camber_pos))) *
               ((1 - 2 * camber_pos) + 2 * camber_pos * x - x * x);
          dyc_dx = (2 * max_camber / ((1 - camber_pos) * (1 - camber_pos))) *
                   (camber_pos - x);
        }
      }

      double theta = std::atan(dyc_dx);

      upper.push_back({x - yt * std::sin(theta), yc + yt * std::cos(theta)});
      lower.push_back({x + yt * std::sin(theta), yc - yt * std::cos(theta)});
    }

    // Closed contour: upper side from leading to trailing edge, then lower
    // side from trailing back to leading edge.
    cached_polygon = upper;
    for (auto it = lower.rbegin(); it != lower.rend(); ++it)
      cached_polygon.push_back(*it);
  }

  // Point-in-polygon test (ray casting) in normalised coordinates.
  bool containsNormalized(double xn, double yn) const {
    bool inside = false;
    size_t n = cached_polygon.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
      double xi = cached_polygon[i].first, yi = cached_polygon[i].second;
      double xj = cached_polygon[j].first, yj = cached_polygon[j].second;
      if (((yi > yn) != (yj > yn)) &&
          (xn < (xj - xi) * (yn - yi) / (yj - yi) + xi)) {
        inside = !inside;
      }
    }
    return inside;
  }

public:
  Airfoil(const types::Coordinate<dim> position_, const double chord_,
          const double thickness_, const double max_camber_ = 0.0,
          const double camber_pos_ = 0.0, const double aoa_deg_ = 0.0)
      : position(position_), chord(chord_), thickness(thickness_),
        max_camber(max_camber_), camber_pos(camber_pos_),
        aoa_rad(aoa_deg_ * 3.14159265358979323846 / 180.0) {
    buildPolygon();
  }

  bool contains(const types::Coordinate<dim> &point) const {
    // trasforma point in coordinate normalizzate locali (inversa di toGrid)
    types::Coordinate<dim> rel = point - position;
    double xr = static_cast<double>(rel.x);
    double yr = static_cast<double>(rel.y);

    // rotazione inversa corretta: stessa forma di toGrid()
    double xs = xr * std::cos(aoa_rad) + yr * std::sin(aoa_rad);
    double ys = -xr * std::sin(aoa_rad) + yr * std::cos(aoa_rad);
    double xn = xs / chord;
    double yn = ys / chord;
    return containsNormalized(xn, yn);
  }

  /// The same box the old getPerimeter() used to scan: one chord per side,
  /// plus two cells of margin for rotation and camber.
  AABB<dim> aabb() const {
    const int half = static_cast<int>(std::ceil(chord)) + 2;
    return {{position.x - half, position.y - half},
            {position.x + half, position.y + half}};
  }
};

} // namespace CollisionDetection
} // namespace lbm

#endif // __SHAPE_HPP
