#pragma once

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/types/common.hpp"

#include <cstddef>

namespace lbm::detail {

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline types::Direction<dim> direction(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_dir<dim, VelocitySet>[diridx];
#else
  return VelocitySet::dir[diridx];
#endif
}

template <typename VelocitySet>
LBM_HD_FUNC inline double weight(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_wi<VelocitySet>[diridx];
#else
  return VelocitySet::wi[diridx];
#endif
}

template <typename VelocitySet>
LBM_HD_FUNC inline std::size_t opposite(const std::size_t diridx) {
#if defined(__CUDA_ARCH__)
  return cuda::vs_opp<VelocitySet>[diridx];
#else
  return VelocitySet::opp[diridx];
#endif
}

} // namespace lbm::detail
