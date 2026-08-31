#pragma once

#include "lbm-sim/backend/cuda/utils.cuh"
#include "lbm-sim/backend/fwd.hpp"
#include "lbm-sim/types/base.hpp"

#include <cuda_runtime.h>

// C++ STD LIB
#include <sstream>
#include <stdexcept>

namespace lbm::profiling {

/// CUDA knobs. Same contract as `BackendProperties<OPEN_MP>`: fill it once at
/// start-up, read it afterwards. The accessors do no locking, so mutating them
/// while a solver runs is a data race.
///
/// Unlike the OpenMP knobs, these do not all live in the runtime. They split
/// in two:
///
///   * *context state* -- the device and its scheduling flags. `apply()`
///     pushes these into the CUDA runtime and `Scope` puts back what was
///     there, exactly like the OpenMP class does with its ICVs.
///
///   * *launch parameters* -- the block shape and the benchmark switch. The
///     runtime has nowhere to store these; the solver reads them straight off
///     this object at launch time. So they take effect as soon as they are
///     set, `Scope` or no `Scope`, and leaving a `scopedApply()` does not undo
///     them. Call `clearBlock()` if you need the default shape back.
template <> class BackendProperties<CUDA> {
public:
  /// The single instance. Same reasoning as the OpenMP specialization: a
  /// function-local static is thread-safely initialized on first use and is
  /// the same object in every translation unit.
  static BackendProperties &get() {
    static BackendProperties instance;
    return instance;
  }

  // A singleton must not be duplicated: a copy would be an object whose
  // mutations nobody else can observe.
  BackendProperties(const BackendProperties &) = delete;
  BackendProperties &operator=(const BackendProperties &) = delete;
  BackendProperties(BackendProperties &&) = delete;
  BackendProperties &operator=(BackendProperties &&) = delete;

  int getDevice() const noexcept { return device_id; }
  unsigned int getSyncPolicy() const noexcept { return sync_policy; }
  bool getBenchmarkMode() const noexcept { return benchmark_mode; }

  /// The block shape to launch with, validated against the selected device.
  ///
  /// Left unset it yields what the solver hardcoded before this class existed
  /// -- 16x16 in 2D, 8x8x8 in 3D -- so an unconfigured process launches
  /// exactly the same grid as it used to.
  ///
  /// A note on picking one: `Grid::field_index` is `ndir * (nx*y + x) + dir`,
  /// an array-of-structures layout, so consecutive threads along x read
  /// `ndir` doubles apart rather than back to back. Coalescing is imperfect
  /// whatever the block shape, but a warp that covers 32 consecutive x still
  /// touches one contiguous span, while 16x16 splits every warp across two
  /// rows and so across two disjoint spans. That makes `block.x` a multiple
  /// of the warp size the first thing worth sweeping.
  template <types::dim_t dim> dim3 getBlock() const {
    const dim3 b =
        isBlockSet() ? block : cuda::create_block_of<dim>(default_extent(dim));
    validateBlock<dim>(b);
    return b;
  }

  /// Which GPU to run on. Only honoured through `apply()`/`scopedApply()`,
  /// which must wrap *all* the device work -- the solver allocates its
  /// `DeviceBuffer`s inside `solve()`, and an allocation lands on whatever
  /// device is current when it runs.
  void setDevice(const int device_id_) noexcept {
    device_id = device_id_;
    // The cached `cudaDeviceProp` describes the old device.
    prop_cache_device = -1;
  }

  void setBlock(const dim3 block_) noexcept { block = block_; }

  void setBlock(const unsigned int x, const unsigned int y,
                const unsigned int z = 1) noexcept {
    block = dim3(x, y, z);
  }

  /// Back to the per-dimension default shape.
  void clearBlock() noexcept { block = dim3(0, 0, 0); }

  /// How the host thread waits on the device: `cudaDeviceScheduleSpin`,
  /// `cudaDeviceScheduleYield`, `cudaDeviceScheduleBlockingSync`, or
  /// `cudaDeviceScheduleAuto` (the default, and what the runtime would have
  /// chosen on its own).
  ///
  /// It shows up in a timing run through every `cudaStreamSynchronize` and
  /// `cudaEventSynchronize`: spinning shaves the wake-up latency off each of
  /// them at the cost of a busy core, blocking gives the core back and pays
  /// the latency. Auto decides by comparing the number of contexts against
  /// the number of logical processors, which is to say the choice can change
  /// under you when something else on the machine starts using the GPU --
  /// pin it if runs have to be comparable.
  void setSyncPolicy(const unsigned int sync_policy_) noexcept {
    sync_policy = sync_policy_ & cudaDeviceScheduleMask;
  }

  /// Turns the iteration loop into a pure compute loop: no periodic
  /// `cudaStreamSynchronize`, no download of rho/u, no frame written, for any
  /// iteration but the last.
  ///
  /// Those three are host-sync and I/O costs, not solver throughput, and
  /// leaving them inside the timed region is what makes a MLUPS number
  /// meaningless. The final iteration still computes the macroscopic fields,
  /// and `solve()` still downloads them once the loop is over, so the lattice
  /// a benchmarked run leaves behind is the same one a normal run would --
  /// only the intermediate frames are gone.
  void setBenchmarkMode(const bool benchmark_mode_) noexcept {
    benchmark_mode = benchmark_mode_;
  }

  /// Cached `cudaGetDeviceProperties` for the selected device. Invalidated by
  /// `setDevice()`.
  const cudaDeviceProp &deviceProperties() const {
    if (prop_cache_device != device_id) {
      LBM_CUDA_CHECK(cudaGetDeviceProperties(&prop_cache, device_id));
      prop_cache_device = device_id;
    }
    return prop_cache;
  }

  /// Push the context state into the CUDA runtime.
  ///
  /// Order matters: `cudaSetDeviceFlags` records the flags of whichever
  /// device is *current*, and initializes it if it is not initialized yet, so
  /// selecting the device has to come first.
  ///
  /// Not `noexcept`, unlike the OpenMP counterpart: asking for a device that
  /// is not there is a real error and `LBM_CUDA_CHECK` throws on it.
  void apply() const {
    LBM_CUDA_CHECK(cudaSetDevice(device_id));
    LBM_CUDA_CHECK(cudaSetDeviceFlags(sync_policy));
  }

  /// RAII counterpart of `apply()`: snapshots the current device and its
  /// scheduling flags, applies the stored ones, and restores the snapshot on
  /// scope exit.
  ///
  ///     {
  ///       const auto scope = BackendProperties<CUDA>::get().scopedApply();
  ///       simulation.solve(solver);
  ///     } // current device and its flags are back to what they were
  ///
  /// Only the context state is restored -- see the class comment.
  class [[nodiscard]] Scope {
  public:
    explicit Scope(const BackendProperties &props) : prev_device(0) {
      LBM_CUDA_CHECK(cudaGetDevice(&prev_device));
      LBM_CUDA_CHECK(cudaSetDevice(props.getDevice()));

      // Read the flags *after* switching: they belong to the device we are
      // about to overwrite them on, which is the target, not `prev_device`.
      unsigned int flags = 0;
      LBM_CUDA_CHECK(cudaGetDeviceFlags(&flags));
      // cudaGetDeviceFlags can hand back cudaDeviceMapHost, which
      // cudaSetDeviceFlags refuses. Keep only the schedule bits -- the only
      // ones we touch.
      prev_flags = flags & cudaDeviceScheduleMask;

      LBM_CUDA_CHECK(cudaSetDeviceFlags(props.getSyncPolicy()));
    }

    ~Scope() {
      // Still current on the target device, so this restores *its* flags;
      // then hand the thread back to the device it came from.
      cudaSetDeviceFlags(prev_flags);
      cudaSetDevice(prev_device);
    }

    // Restoring twice, or from the wrong scope, would defeat the purpose.
    Scope(const Scope &) = delete;
    Scope &operator=(const Scope &) = delete;
    Scope(Scope &&) = delete;
    Scope &operator=(Scope &&) = delete;

  private:
    int prev_device;
    unsigned int prev_flags;
  };

  /// Returned by value: C++17 guaranteed copy elision builds the `Scope`
  /// straight into the caller's variable, so the deleted move is never
  /// needed.
  [[nodiscard]] Scope scopedApply() const { return Scope(*this); }

private:
  /// Defaults reproduce the behaviour the solver had before this class
  /// existed: device 0, the block shape that was hardcoded at the launch
  /// site, the runtime's own scheduling heuristic, and frames written as
  /// usual.
  BackendProperties() noexcept
      : device_id(0), block(0, 0, 0), sync_policy(cudaDeviceScheduleAuto),
        benchmark_mode(false), prop_cache(), prop_cache_device(-1) {}

  ~BackendProperties() = default;

  /// `dim3` has no "unset" value, and a block with a zero extent could never
  /// be launched anyway, so all-zero doubles as the sentinel.
  bool isBlockSet() const noexcept { return block.x != 0; }

  static unsigned int default_extent(const types::dim_t dim) noexcept {
    return dim == 2 ? 16 : 8;
  }

  /// Structural checks first, device limits second: the former need no
  /// driver, so a block that is wrong on its own terms is reported as such
  /// even where there is no GPU to query.
  template <types::dim_t dim> void validateBlock(const dim3 &b) const {
    // In 2D the kernels index with `thread_coordinate<2>()`, which ignores
    // threadIdx.z: a block deeper than one plane would hand the same cell to
    // several threads, and they would all write it. That is a race, not a
    // slow launch, so it has to be caught here rather than left to the
    // runtime -- which would happily accept it.
    if constexpr (dim == 2) {
      if (b.z != 1) {
        throw std::invalid_argument(
            "BackendProperties<CUDA>: a 2D block must have z == 1; a deeper "
            "block maps several threads onto the same cell.");
      }
    }

    if (b.x == 0 || b.y == 0 || b.z == 0) {
      throw std::invalid_argument(
          "BackendProperties<CUDA>: no block extent may be zero.");
    }

    const cudaDeviceProp &prop = deviceProperties();

    const unsigned int threads = b.x * b.y * b.z;
    if (threads > static_cast<unsigned int>(prop.maxThreadsPerBlock)) {
      std::ostringstream oss;
      oss << "BackendProperties<CUDA>: block " << b.x << "x" << b.y << "x"
          << b.z << " asks for " << threads << " threads, but device "
          << device_id << " (" << prop.name << ") allows at most "
          << prop.maxThreadsPerBlock << " per block.";
      throw std::invalid_argument(oss.str());
    }

    const unsigned int extent[3] = {b.x, b.y, b.z};
    for (int d = 0; d < 3; ++d) {
      if (extent[d] > static_cast<unsigned int>(prop.maxThreadsDim[d])) {
        std::ostringstream oss;
        oss << "BackendProperties<CUDA>: block extent " << extent[d]
            << " along axis " << d << " exceeds the " << prop.maxThreadsDim[d]
            << " allowed by device " << device_id << " (" << prop.name << ").";
        throw std::invalid_argument(oss.str());
      }
    }
  }

  int device_id;
  dim3 block;
  unsigned int sync_policy;
  bool benchmark_mode;

  mutable cudaDeviceProp prop_cache;
  mutable int prop_cache_device;
};

} // namespace lbm::profiling
