#ifndef __LBM_SIM_BACKEND_HPP
#define __LBM_SIM_BACKEND_HPP

#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/cuda/annotations.hpp"

namespace lbm {

enum ExecutionBackend { CUDA, OPEN_MP };

namespace detail { // using cuda or openmp

template <unsigned short int dim, typename VelocitySet>
LBM_HD_FUNC inline types::VectorInt<dim> direction(const std::size_t diridx) {
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

} // namespace detail

} // namespace lbm

#endif // __LBM_SIM_BACKEND_HPP
