#ifndef __LBM_SIM_BACKEND_METADATA_HPP
#define __LBM_SIM_BACKEND_METADATA_HPP

namespace lbm {

enum ExecutionBackend { CUDA, OPEN_MP };

template <ExecutionBackend backend> struct ExecutionContext;

template <> struct ExecutionContext<ExecutionBackend::OPEN_MP> {};

template <> struct ExecutionContext<ExecutionBackend::CUDA> {
  // Opaque stream handle to avoid forcing CUDA headers in common code.
  void *stream = nullptr;
  int block_x = 16;
  int block_y = 16;
  int block_z = 1;
};

} // namespace lbm

#endif // __LBM_SIM_BACKEND_METADATA_HPP
