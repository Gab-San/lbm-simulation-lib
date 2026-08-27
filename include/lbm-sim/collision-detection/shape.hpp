#ifndef __SHAPE_HPP
#define __SHAPE_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/collision-detection/algorithms/collision.hpp"
#include "lbm-sim/collision-detection/algorithms/rasterization.hpp"

// C++ STANDARD LIB
#include <stdexcept>
#include <vector>

// FIXME: Should be extendible to 3D

namespace lbm {
namespace CollisionDetection {

template <types::dim_t dim, class DerivedShape> class Shape {
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

template <types::dim_t dim> class Segment : public Shape<dim, Segment<dim>> {
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

template <types::dim_t dim> class Circle : public Shape<dim, Circle<dim>> {
  const types::Coordinate<dim> center;
  const unsigned int radius;
  mutable std::vector<types::Coordinate<dim>> cached_perimeter;

public:
  Circle(const types::Coordinate<dim> center_, const unsigned int radius_)
      : center(center_), radius(radius_) {}

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    if constexpr (dim == 2) {
      // return checkCollisionWith(point);
      return contains(point);
    } else {
      return checkCollisionWith(point);
    }
  }

  bool contains(const types::Coordinate<dim> &point) const {
    const int center_point_dist_x = center.x - point.x;
    const unsigned int dist_x = center_point_dist_x * center_point_dist_x;
    const int center_point_dist_y = center.y - point.y;
    const unsigned int dist_y = center_point_dist_y * center_point_dist_y;

    // Contains so we use "<=" for the inner points
    return dist_x + dist_y <= radius * radius;
  }

  // std::vector<types::Coordinate<dim>> &getPerimeter() const {
  //   if (!cached_perimeter.empty()) {
  //     return cached_perimeter;
  //   }
  //   cached_perimeter.clear();
  //   int x = -radius;
  //   int y = 0;
  //   cached_perimeter.push_back(center + types::Coordinate<dim>(x, y));
  //   cached_perimeter.push_back(center + types::Coordinate<dim>(-x, y));
  //   cached_perimeter.push_back(center + types::Coordinate<dim>(0, radius));
  //   cached_perimeter.push_back(center + types::Coordinate<dim>(0, -radius));
  //   while (x < static_cast<int>(radius)) {
  //     if (x < 0) {
  //       if (x * x + (y + 1) * (y + 1) > radius * radius) {
  //         x++;
  //       } else {
  //         y++;
  //       }
  //     } else {
  //       if ((x + 1) * (x + 1) + y * y > radius * radius) {
  //         y--;
  //       } else {
  //         x++;
  //       }
  //     }
  //     if (y == -y) {
  //       continue;
  //     }
  //     cached_perimeter.push_back(center + types::Coordinate<dim>(x, y));
  //     cached_perimeter.push_back(center + types::Coordinate<dim>(x, -y));
  //   }
  //   return cached_perimeter;
  // }

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    if (!cached_perimeter.empty()) {
      return cached_perimeter;
    }

    cached_perimeter.clear();

    const int r = static_cast<int>(radius);
    const unsigned int r2 = radius * radius;

    for (int dx = -r; dx <= r; ++dx) {
      for (int dy = -r; dy <= r; ++dy) {
        if (static_cast<unsigned int>(dx * dx + dy * dy) <= r2) {
          cached_perimeter.push_back(center + types::Coordinate<dim>(dx, dy));
        }
      }
    }

    return cached_perimeter;
  }

private:
  bool checkCollisionWith(const types::Coordinate<2> &point) const {
    const int center_point_dist_x = center.x - point.x;
    const unsigned int dist_x = center_point_dist_x * center_point_dist_x;
    const int center_point_dist_y = center.y - point.y;
    const unsigned int dist_y = center_point_dist_y * center_point_dist_y;

    return dist_x + dist_y <= radius * radius;
  }

  bool checkCollisionWith(const types::Coordinate<3> &point) const {
    throw std::runtime_error("Circle::checkCollisionWith(utils::Point<3>&) : "
                             "Collision Detection in 3D not yet implemented!");
  }
};

template <types::dim_t dim>
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

template <unsigned short int dim>
class Airfoil : public Shape<dim, Airfoil<dim>> {
  static_assert(dim == 2, "Airfoil currently supports only 2D");

  const types::Coordinate<dim> position; // punto di riferimento (leading edge)
  const double chord;                    // corda in celle
  const double thickness;                // XX/100 (es. 0.12 per NACA 0012)
  const double max_camber;               // M/100 (es. 0.02 per NACA 2412)
  const double camber_pos;               // P/10 (es. 0.4 per NACA 2412)
  const double aoa_rad;                  // angolo di attacco in radianti

  mutable std::vector<types::Coordinate<dim>> cached_perimeter;
  mutable std::vector<std::pair<double, double>>
      cached_polygon; // punti normalizzati (x,y) del contorno chiuso

  // genera il contorno chiuso (upper poi lower al contrario) in coordinate
  // normalizzate [0,1] di corda
  void buildPolygon() const {
    if (!cached_polygon.empty())
      return;

    constexpr int N = 100; // risoluzione: piu' alto = piu' preciso
    std::vector<std::pair<double, double>> upper, lower;

    for (int i = 0; i <= N; ++i) {
      // cosine spacing per infittire vicino a leading/trailing edge
      // definisco la costante M_PI per sicurezza, anche se dovrebbe essere
      // definita in cmath
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

    // contorno chiuso: upper da leading a trailing, poi lower da trailing a
    // leading
    cached_polygon = upper;
    for (auto it = lower.rbegin(); it != lower.rend(); ++it)
      cached_polygon.push_back(*it);
  }

  // ruota e scala un punto normalizzato -> coordinate assolute di griglia
  types::Coordinate<dim> toGrid(double xn, double yn) const {
    // scala per la corda
    double xs = xn * chord;
    double ys = yn * chord;
    // ruota per angolo di attacco (attorno al leading edge, origine locale)
    double xr = xs * std::cos(aoa_rad) + ys * std::sin(aoa_rad);
    double yr = -xs * std::sin(aoa_rad) + ys * std::cos(aoa_rad);
    return position + types::Coordinate<dim>(static_cast<int>(std::round(xr)),
                                             static_cast<int>(std::round(yr)));
  }

  // point-in-polygon (ray casting) in coordinate normalizzate
  bool containsNormalized(double xn, double yn) const {
    buildPolygon();
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
        aoa_rad(aoa_deg_ * 3.14159265358979323846 / 180.0) {}

  bool contains(const types::Coordinate<dim> &point) const {
    // trasforma point in coordinate normalizzate locali (inversa di toGrid)
    types::Coordinate<dim> rel = point - position;
    double xr = static_cast<double>(rel.x);
    double yr = static_cast<double>(rel.y);
    // rotazione inversa
    // double xs = xr * std::cos(aoa_rad) - yr * std::sin(aoa_rad);
    // double ys = xr * std::sin(aoa_rad) + yr * std::cos(aoa_rad);

    // rotazione inversa corretta: stessa forma di toGrid()
    double xs = xr * std::cos(aoa_rad) + yr * std::sin(aoa_rad);
    double ys = -xr * std::sin(aoa_rad) + yr * std::cos(aoa_rad);
    double xn = xs / chord;
    double yn = ys / chord;
    return containsNormalized(xn, yn);
  }

  bool isCollidingWith(const types::Coordinate<dim> &point) const {
    return contains(point);
  }

  std::vector<types::Coordinate<dim>> &getPerimeter() const {
    if (!cached_perimeter.empty())
      return cached_perimeter;
    buildPolygon();

    cached_perimeter.clear();

    // bounding box in coordinate di griglia (con margine per rotazione/camber)
    int half = static_cast<int>(std::ceil(chord)) + 2;
    for (int dx = -half; dx <= half; ++dx) {
      for (int dy = -half; dy <= half; ++dy) {
        const types::Coordinate<dim> p =
            position + types::Coordinate<dim>(dx, dy);
        if (contains(p)) {
          cached_perimeter.push_back(p);
        }
      }
    }

    return cached_perimeter;
  }
};

} // namespace CollisionDetection
} // namespace lbm

#endif // __SHAPE_HPP
