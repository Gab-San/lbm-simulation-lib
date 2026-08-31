#ifndef __LBM_SIM_CUDA_STRUCTS_CUH
#define __LBM_SIM_CUDA_STRUCTS_CUH

#include "lbm-sim/cuda/utils.cuh"

#include <cuda_runtime.h>

#include <assert.h>
#include <vector>

namespace lbm {
namespace cuda {
template <typename T> class DeviceBuffer {
  T *_ptr;
  std::size_t _size;

public:
  DeviceBuffer() : _ptr(nullptr), _size(0) {}

  explicit DeviceBuffer(std::size_t cap_) : _size(cap_) {
    LBM_CUDA_CHECK(cudaMalloc(&_ptr, _size * sizeof(T)));
  }

  ~DeviceBuffer() {
    if (_ptr)
      cudaFree(_ptr);
  }

  DeviceBuffer(DeviceBuffer &&o) noexcept
      : _ptr(std::exchange(o._ptr, nullptr)), _size(std::exchange(o._size, 0)) {
  }

  DeviceBuffer &operator=(DeviceBuffer &&o) noexcept {
    if (this != &o) {
      if (_ptr)
        cudaFree(_ptr);
      _ptr = std::exchange(o._ptr, nullptr);
      _size = std::exchange(o._size, 0);
    }
    return *this;
  }

  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer &operator=(const DeviceBuffer &) = delete;

  void swap(DeviceBuffer &o) noexcept {
    std::swap(_ptr, o._ptr);
    std::swap(_size, o._size);
  }

  T *data() const noexcept { return _ptr; }
  std::size_t size() const noexcept { return _size; }
  std::size_t bytes() const noexcept { return _size * sizeof(T); }

  void upload_async(const T *src, std::size_t cpy_size, cudaStream_t stream) {
    assert(cpy_size <= _size && "DeviceBuffer<T>::upload_async(): upload would "
                                "overflow device buffer.");
    LBM_CUDA_CHECK(cudaMemcpyAsync(this->_ptr, src, cpy_size * sizeof(T),
                                   cudaMemcpyHostToDevice, stream));
  }

  template <typename Container>
  void upload_async(const Container &src, cudaStream_t stream) {
    upload_async(src.data(), src.size(), stream);
  }

  void download_async(T *dst, std::size_t cpy_size, cudaStream_t stream) const {
    assert(cpy_size <= _size &&
           "DeviceBuffer<T>::download_async(): download would "
           "overflow host buffer.");
    LBM_CUDA_CHECK(cudaMemcpyAsync(dst, this->_ptr, cpy_size * sizeof(T),
                                   cudaMemcpyDeviceToHost, stream));
  }

  template <typename Container>
  void download_async(Container &dst, cudaStream_t stream) const {
    download_async(dst.data(), dst.size(), stream);
  }
};

class StreamHandler {
private:
  cudaStream_t s;

public:
  StreamHandler() : s(nullptr) {}

  ~StreamHandler() {
    if (s)
      cudaStreamDestroy(s);
  }

  StreamHandler(const StreamHandler &) = delete;
  StreamHandler &operator=(const StreamHandler &) = delete;

  operator cudaStream_t() const noexcept {
    return s;
  } // implicit, so it drops into launches
};
} // namespace cuda
} // namespace lbm

#endif // __LBM_SIM_CUDA_STRUCTS_CUH
