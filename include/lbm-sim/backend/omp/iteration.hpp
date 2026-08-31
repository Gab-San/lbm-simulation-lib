#ifndef __LBM_SIM_OMP_ITERATION_HPP
#define __LBM_SIM_OMP_ITERATION_HPP

#include "lbm-sim/backend/omp/annotations.hpp"
#include "lbm-sim/types/common.hpp"

namespace lbm {
namespace iteration {
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
