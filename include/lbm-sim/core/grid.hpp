/**
 * @file grid.hpp
 * @brief Grid<dim>: the domain extents plus the index arithmetic every
 *        backend agrees on.
 *
 * A Grid carries no data. It is a size and the four functions that turn a
 * coordinate into an offset, small enough to be passed by value into a CUDA
 * kernel and cheap enough to be recomputed rather than cached.
 *
 * Two layouts are defined here and used everywhere else:
 * - **scalar fields** (@c rho, @c u, the velocity norms) are indexed by
 *   scalar_index(), with x varying fastest;
 * - **population fields** are indexed by field_index(), which is
 *   direction-major *within* a cell: all @c ndir populations of one node sit
 *   contiguously, so a collision step touches one cache line per node.
 */

#pragma once

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/backend/cuda/annotations.hpp"

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {

/**
 * @brief Domain extents and the index arithmetic over them.
 *
 * @tparam dim Spatial dimension (2 or 3).
 *
 * @note Trivially copyable and @c const-only, so it can be captured by value
 *       in a device lambda or passed straight as a kernel argument.
 */
template <types::dim_t dim> struct Grid {
  /// Number of cells along each axis. @c z is absent when @c dim == 2.
  const types::DimPoint<dim> size;

  /**
   * @brief Builds a grid of the given extents.
   * @param grid_dim_ Cell counts per axis.
   */
  LBM_HD_FUNC Grid(types::DimPoint<dim> grid_dim_) : size(grid_dim_) {}

  /**
   * @brief Offset of a cell in a scalar field, x varying fastest.
   *
   * @f$ i = x + N_x (y + N_y z) @f$, with the @c z term dropped in 2D.
   *
   * @param p Cell coordinate.
   * @return Index into a field of getArea() entries.
   *
   * @warning No bounds check. Call contains() first when @p p may be outside
   *          the domain; boundary handling relies on that being the caller's
   *          job.
   */
  LBM_HD_FUNC inline std::size_t
  scalar_index(const types::Coordinate<dim> &p) const {
    if constexpr (dim == 2) {
      return size.x * p.y + p.x;
    } else {
      return size.x * (size.y * p.z + p.y) + p.x;
    }
  }

  /**
   * @brief Offset of one population of one cell.
   *
   * Equivalent to @c ndir*scalar_index(p) + dir: the @p ndir populations of
   * a node are contiguous, and @p dir is the offset inside that group.
   *
   * @param p    Cell coordinate.
   * @param dir  Direction index, in @c [0, ndir).
   * @param ndir Number of directions of the velocity set in use.
   * @return Index into a field of @c getArea()*ndir entries.
   *
   * @note @p ndir is a parameter rather than a template argument so that the
   *       grid stays independent of the velocity set; callers pass
   *       @c VelocitySet::ndir, which is @c constexpr, so the multiply is
   *       folded away.
   */
  LBM_HD_FUNC inline std::size_t field_index(const types::Coordinate<dim> &p,
                                             std::size_t dir,
                                             std::size_t ndir) const {
    return ndir * scalar_index(p) + dir;
  }

  /**
   * @brief Total cell count: @c nx*ny in 2D, @c nx*ny*nz in 3D.
   *
   * The name is historical -- it is a volume in 3D -- and is kept because it
   * appears in every field allocation in the library.
   */
  LBM_HD_FUNC inline std::size_t getArea() const {
    if constexpr (dim == 2) {
      return size.x * size.y;
    } else {
      return size.x * size.y * size.z;
    }
  }

  /**
   * @brief The extents as a plain array, for the generic iteration helpers.
   * @see lbm::iteration::unflatten()
   */
  inline std::array<std::size_t, dim> extents() const {
    if constexpr (dim == 2) {
      return {size.x, size.y};
    } else {
      return {size.x, size.y, size.z};
    }
  }

  /**
   * @brief True when @p p lies inside the domain on every axis.
   *
   * The bound is compared as @c int on purpose: coordinates are signed, and
   * link resolution routinely produces a negative component before deciding
   * what the boundary does with it.
   */
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
