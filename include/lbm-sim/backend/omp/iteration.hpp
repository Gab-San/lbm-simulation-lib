/**
 * @file iteration.hpp
 * @brief Turning a flat cell index back into a coordinate.
 *
 * The host solver parallelises over a single flattened loop rather than a
 * nested one: one @c omp @c parallel @c for over @c getArea() iterations
 * gives the runtime one chunkable range and behaves identically in 2D and
 * 3D. Each iteration then needs its coordinate back, which is what
 * unflatten() is for.
 */

#ifndef __LBM_SIM_OMP_ITERATION_HPP
#define __LBM_SIM_OMP_ITERATION_HPP

#include "lbm-sim/backend/omp/annotations.hpp"
#include "lbm-sim/types/common.hpp"

namespace lbm {
namespace iteration {

/**
 * @brief Inverse of Grid::scalar_index(): flat index to coordinate.
 *
 * @param idx Flat index, in @c [0, ext[0]*ext[1]*...).
 * @param ext Domain extents, from Grid::extents().
 * @return The coordinate that scalar_index() maps back to @p idx.
 *
 * x varies fastest, matching the scalar field layout, so consecutive
 * iterations of the parallel loop touch consecutive memory.
 *
 * @warning No bounds check: an @p idx past the end silently produces an
 *          out-of-domain coordinate.
 */
template <types::dim_t dim>
inline types::Coordinate<dim>
unflatten(std::size_t idx, const std::array<std::size_t, dim> &ext) noexcept {
  std::array<int, dim> c;

  UNROLL_FULL
  for (unsigned short d = 0; d < dim; ++d) {
    c[d] = static_cast<int>(idx % static_cast<std::size_t>(ext[d]));
    idx /= static_cast<std::size_t>(ext[d]);
  }

  if constexpr (dim == 2)
    return types::Coordinate<2>(c[0], c[1]);
  else
    return types::Coordinate<3>(c[0], c[1], c[2]);
}
} // namespace iteration
} // namespace lbm

#endif // __LBM_SIM_OMP_ITERATION_HPP
