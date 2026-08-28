#ifndef __LBM_SIM_BACKEND_CUDA_UTILS_CUH
#define __LBM_SIM_BACKEND_CUDA_UTILS_CUH

#include "lbm-sim/types/base.hpp"
#include "lbm/logging.hpp"

#include "quill/LogMacros.h"

#include <cuda_runtime.h>

// C++ STD LIB
#include <sstream>
#include <stdexcept>

namespace lbm {

// Error checks for async errors
#define LBM_CUDA_CHECK(expr)                                                   \
  do {                                                                         \
    const cudaError_t _lbm_err = (expr);                                       \
    if (_lbm_err != cudaSuccess) {                                             \
      std::ostringstream _lbm_oss;                                             \
      _lbm_oss << "CUDA ERROR " << cudaGetErrorString(_lbm_err) << " at "      \
               << __FILE__ << " : " << __LINE__;                               \
      throw std::runtime_error(_lbm_oss.str());                                \
    }                                                                          \
  } while (0)

namespace cuda {

template <types::dim_t dim>
__device__ inline types::Coordinate<dim> thread_coordinate() {
  if constexpr (dim == 2) {
    return types::Coordinate<2>(blockIdx.x * blockDim.x + threadIdx.x,
                                blockIdx.y * blockDim.y + threadIdx.y);
  } else {
    return types::Coordinate<3>(blockIdx.x * blockDim.x + threadIdx.x,
                                blockIdx.y * blockDim.y + threadIdx.y,
                                blockIdx.z * blockDim.z + threadIdx.z);
  }
}

inline unsigned int ceil_div(unsigned int a, unsigned int b) {
  return (a + b - 1) / b;
}

template <types::dim_t dim>
inline dim3 ceil_div(types::DimPoint<dim> size, dim3 b) {
  if constexpr (dim == 2) {
    return dim3(ceil_div(size.x, b.x), ceil_div(size.y, b.y));
  } else {
    return dim3(ceil_div(size.x, b.x), ceil_div(size.y, b.y),
                ceil_div(size.z, b.z));
  }
}

template <types::dim_t dim> inline dim3 create_block_of(unsigned int n) {
  if constexpr (dim == 2) {
    return dim3(n, n);
  } else {
    return dim3(n, n, n);
  }
}

inline void log_device_info(quill::Logger *logger) {
  int deviceCount = 0;
  LBM_CUDA_CHECK(cudaGetDeviceCount(&deviceCount));
  LOG_INFO(logger, "System has {} CUDA device(s).", deviceCount);

  for (int i = 0; i < deviceCount; ++i) {
    cudaDeviceProp prop;
    LBM_CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
    LOG_INFO(logger, "  Device {}: {}", i, prop.name);
    LOG_INFO(logger, "    SMs: {}", prop.multiProcessorCount);
    LOG_INFO(logger, "    Max threads/block: {}", prop.maxThreadsPerBlock);
    LOG_INFO(logger, "    Max threads/SM: {}",
             prop.maxThreadsPerMultiProcessor);
    LOG_INFO(logger, "    Warp size: {}", prop.warpSize);
  }
}

} // namespace cuda

} // namespace lbm

#endif //  __LBM_SIM_BACKEND_CUDA_UTILS_HPP
