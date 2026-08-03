#ifndef __LBM_SIM_CORE_GRID_HPP
#define __LBM_SIM_CORE_GRID_HPP

#include "collision-detection/core/types.hpp"

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {
// constexpr bool always_false = false;
template <unsigned short int dim> struct Grid {
  const CollisionDetection::types::DimPoint<dim> size;

  Grid(CollisionDetection::types::DimPoint<dim> grid_dim_) : size(grid_dim_) {}

  inline const std::size_t
  scalar_index(const CollisionDetection::types::Coordinate<dim> &p) const {
    if constexpr (dim == 2) {
      return size.x * p.y + p.x;
    } else {
      static_assert(dim != dim,
                    "Grid<3>::scalar_index() : 3D not implemented yet!");
      // TODO: Check implementation
      //
      // return Nx * (Ny * z + y) + x;
    }
  }

  // Index position of a cell for a direction defined vector
  // This function is equal to: (Nx*Ny*dir) + (Nx*y)+x
  // making dir work as an offset.
  inline std::size_t
  field_index(const CollisionDetection::types::Coordinate<dim> &p,
              std::size_t dir, std::size_t ndir) const {
    if constexpr (dim == 2) {
      // NOTE: maybe grid could be templated on velocity sets.
      return ndir * (size.x * p.y + p.x) + dir;
    } else {
      static_assert(dim != dim,
                    "Grid<3>::field_index() : 3D not implemented yet!");
      // TODO: Check implementation
      //
      // return Nx * (Ny * (Nz * dir + z) + y) + x;
    }
  }

  inline std::size_t getArea() const {
    if constexpr (dim == 2) {
      return size.x * size.y;
    } else {
      return size.x * size.y * size.z;
    }
  }

  inline bool
  contains(const CollisionDetection::types::Coordinate<dim> &p) const {
    const bool isIn = p.x >= 0 && p.x < size.x && p.y >= 0 && p.y < size.y;
    if constexpr (dim == 2) {
      return isIn;
    } else {
      return isIn && p.z >= 0 && p.z < size.z;
    }
  }
};

} // namespace lbm

#endif // __LBM_SIM_CORE_GRID_HPP
