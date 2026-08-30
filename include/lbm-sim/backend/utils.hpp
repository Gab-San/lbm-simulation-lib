/**
 * @file utils.hpp
 * @brief Velocity-set table access that works on both the host and the
 *        device.
 *
 * A velocity set stores its tables as @c constexpr arrays, which the host
 * reads directly. Device code cannot: the arrays have to live in
 * @c __constant__ memory, uploaded once by
 * lbm::cuda::upload_lattice_constants(). These three accessors hide that
 * split behind a @c __CUDA_ARCH__ switch, so the boundary conditions and the
 * collision kernels are written once against @c direction(), @c weight() and
 * @c opposite() instead of against either storage.
 *
 * @warning On the device the accessors read the @c __constant__ tables
 *          unconditionally. Forgetting upload_lattice_constants() before the
 *          first kernel launch yields zeros, not a diagnostic.
 */

#pragma once

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/types/common.hpp"

#include <cstddef>

namespace lbm::detail {

/// @brief The discrete velocity @f$ \mathbf{c}_i @f$ of direction @p diridx.
template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline types::Direction<dim> direction(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_dir<dim, VelocitySet>[diridx];
#else
  return VelocitySet::dir[diridx];
#endif
}

/// @brief The lattice weight @f$ w_i @f$ of direction @p diridx.
template <typename VelocitySet>
LBM_HD_FUNC inline double weight(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_wi<VelocitySet>[diridx];
#else
  return VelocitySet::wi[diridx];
#endif
}

/// @brief The index of the direction opposite to @p diridx, i.e. the
///        @f$ \bar\imath @f$ with @f$ \mathbf{c}_{\bar\imath} =
///        -\mathbf{c}_i @f$. This is what bounce-back reflects into.
template <typename VelocitySet>
LBM_HD_FUNC inline std::size_t opposite(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_opp<VelocitySet>[diridx];
#else
  return VelocitySet::opp[diridx];
#endif
}

} // namespace lbm::detail
