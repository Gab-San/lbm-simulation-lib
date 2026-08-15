#ifndef __LBM_SIM_BACKEND_CUDA_UTILS_HPP
#define __LBM_SIM_BACKEND_CUDA_UTILS_HPP

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

namespace cuda_detail {

inline unsigned int ceil_div(unsigned int a, unsigned int b) {
  return (a + b - 1) / b;
}

} // namespace cuda_detail

} // namespace lbm

#endif //  __LBM_SIM_BACKEND_CUDA_UTILS_HPP
