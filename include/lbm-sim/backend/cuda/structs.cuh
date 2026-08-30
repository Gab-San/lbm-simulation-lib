/**
 * @file structs.cuh
 * @brief DeviceBuffer: an owning, move-only handle on a device allocation.
 *
 * The device-side counterpart of a @c std::vector, reduced to what the
 * solver needs: allocate, upload, download, free. Move-only, so ownership of
 * the pointer is unambiguous and the deallocation happens exactly once, at
 * scope exit, including on the path where a launch throws.
 */

#ifndef __LBM_SIM_CUDA_STRUCTS_CUH
#define __LBM_SIM_CUDA_STRUCTS_CUH

#include "lbm-sim/backend/cuda/utils.cuh"

#include <cuda_runtime.h>

#include <assert.h>
#include <vector>

namespace lbm {
namespace cuda {
/**
 * @brief An owning, move-only handle on a device allocation.
 *
 * @tparam T Element type; the buffer is sized in elements, not bytes.
 *
 * Copy is deleted and move transfers the pointer, so ownership is
 * unambiguous and @c cudaFree runs exactly once, at scope exit -- including
 * on the path where a launch throws through LBM_CUDA_CHECK.
 *
 * @note swap() is what makes the population double-buffering work: the
 *       solver swaps the two buffers at the end of each step, which is a
 *       pointer exchange and not a copy.
 */
template <typename T> class DeviceBuffer {
  T *_ptr;
  std::size_t _size;

public:
  /// @brief An empty buffer owning nothing.
  DeviceBuffer() : _ptr(nullptr), _size(0) {}

  /// @brief Allocates @p cap_ elements on the current device.
  /// @throws std::runtime_error if the allocation fails.
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

  /// @brief Exchanges the two allocations. No copy, no reallocation.
  void swap(DeviceBuffer &o) noexcept {
    std::swap(_ptr, o._ptr);
    std::swap(_size, o._size);
  }

  /// @brief The raw device pointer, to be passed as a kernel argument.
  /// @warning Not dereferenceable from the host.
  T *data() const noexcept { return _ptr; }

  /// @brief Capacity, in elements.
  std::size_t size() const noexcept { return _size; }

  /// @brief Capacity, in bytes.
  std::size_t bytes() const noexcept { return _size * sizeof(T); }

  /**
   * @brief Host-to-device copy, enqueued on @p stream.
   *
   * @param src      Host source.
   * @param cpy_size Number of elements to copy.
   * @param stream   Stream to enqueue on.
   *
   * @warning Asynchronous: @p src must stay alive and unmodified until the
   *          stream has been synchronised. The overflow check is an
   *          @c assert, so it disappears with @c NDEBUG.
   */
  void upload_async(const T *src, std::size_t cpy_size, cudaStream_t stream) {
    assert(cpy_size <= _size && "DeviceBuffer<T>::upload_async(): upload would "
                                "overflow device buffer.");
    LBM_CUDA_CHECK(cudaMemcpyAsync(this->_ptr, src, cpy_size * sizeof(T),
                                   cudaMemcpyHostToDevice, stream));
  }

  /// @brief Uploads a whole host container (anything with @c data() and
  ///        @c size()).
  template <typename Container>
  void upload_async(const Container &src, cudaStream_t stream) {
    upload_async(src.data(), src.size(), stream);
  }

  /// @brief Device-to-host copy, enqueued on @p stream. Same asynchrony
  ///        caveat as upload_async().
  void download_async(T *dst, std::size_t cpy_size, cudaStream_t stream) const {
    assert(cpy_size <= _size &&
           "DeviceBuffer<T>::download_async(): download would "
           "overflow host buffer.");
    LBM_CUDA_CHECK(cudaMemcpyAsync(dst, this->_ptr, cpy_size * sizeof(T),
                                   cudaMemcpyDeviceToHost, stream));
  }

  /// @brief Downloads into a whole host container.
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
