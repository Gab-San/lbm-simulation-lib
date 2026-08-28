#pragma once

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/cuda/annotations.hpp"

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {

template <types::dim_t dim> struct Grid {
  const types::DimPoint<dim> size;

  LBM_HD_FUNC Grid(types::DimPoint<dim> grid_dim_) : size(grid_dim_) {}

  LBM_HD_FUNC inline std::size_t
  scalar_index(const types::Coordinate<dim> &p) const {
    if constexpr (dim == 2) {
      return size.x * p.y + p.x;
    } else {
      return size.x * (size.y * p.z + p.y) + p.x;
    }
  }

  // Index position of a cell for a direction defined vector
  // This function is equal to: (Nx*Ny*dir) + (Nx*y)+x
  // making dir work as an offset.
  LBM_HD_FUNC inline std::size_t field_index(const types::Coordinate<dim> &p,
                                             std::size_t dir,
                                             std::size_t ndir) const {
    return ndir * scalar_index(p) + dir;
  }

  LBM_HD_FUNC inline std::size_t getArea() const {
    if constexpr (dim == 2) {
      return size.x * size.y;
    } else {
      return size.x * size.y * size.z;
    }
  }

  inline std::array<std::size_t, dim> extents() const {
    if constexpr (dim == 2) {
      return {size.x, size.y};
    } else {
      return {size.x, size.y, size.z};
    }
  }

  LBM_HD_FUNC inline bool contains(const types::Coordinate<dim> &p) const {
    const bool isIn = p.x >= 0 && p.x < static_cast<int>(size.x) && p.y >= 0 &&
                      p.y < static_cast<int>(size.y);
    if constexpr (dim == 2)
      return isIn;
    else
      return isIn && p.z >= 0 && p.z < static_cast<int>(size.z);
  }
};

} // namespace lbm
