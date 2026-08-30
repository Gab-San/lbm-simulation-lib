/**
 * @file utils.cuh
 * @brief CUDA error checking and the stream/launch helpers.
 *
 * @c LBM_CUDA_CHECK() wraps a runtime call and turns a non-success status
 * into a @c std::runtime_error naming the call site. It is used on every
 * allocation, copy and launch: CUDA errors are sticky and otherwise surface
 * far from their cause.
 */

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

/// @brief Runs a CUDA runtime call and turns a failure into a
///        @c std::runtime_error naming the file and line.
///
/// Used on every allocation, copy and launch: CUDA errors are sticky, so an
/// unchecked failure surfaces later and somewhere else.
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

/// @brief The lattice node this thread is responsible for, from its block
///        and thread indices.
/// @warning The launch grid is rounded up, so the result may lie outside the
///          domain: every kernel checks Grid::contains() before doing work.
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

/// @brief Integer division rounding up.
inline unsigned int ceil_div(unsigned int a, unsigned int b) {
  return (a + b - 1) / b;
}

/// @brief The launch grid that covers @p size with blocks of shape @p b,
///        rounding up on every axis.
template <types::dim_t dim>
inline dim3 ceil_div(types::DimPoint<dim> size, dim3 b) {
  if constexpr (dim == 2) {
    return dim3(ceil_div(size.x, b.x), ceil_div(size.y, b.y));
  } else {
    return dim3(ceil_div(size.x, b.x), ceil_div(size.y, b.y),
                ceil_div(size.z, b.z));
  }
}

/// @brief A cubic block of side @p n: @c n*n threads in 2D, @c n*n*n in 3D.
/// @warning The product must stay within the device's maximum threads per
///          block -- 1024 on every architecture the project targets, so
///          @c n <= 10 in 3D.
template <types::dim_t dim> inline dim3 create_block_of(unsigned int n) {
  if constexpr (dim == 2) {
    return dim3(n, n);
  } else {
    return dim3(n, n, n);
  }
}

/// @brief Logs the visible devices and their headline limits, once per run.
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
