#ifndef __LBM_SIM_SOLVER_SOLVER_CUDA_CUH
#define __LBM_SIM_SOLVER_SOLVER_CUDA_CUH

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/backend/cuda-utils.cuh"
#include "lbm-sim/backend/metadata.hpp"

// CUDA LIB
#include <cuda_runtime.h>

// C++ STD LIB
#include <utility>

// WARN: LBM_CUDA_CHECK is in backend/cuda-utils.cuh

namespace lbm {
namespace cuda_detail {}

template <enum CollisionModel cm_t>
class TMPCUDASolver2D : public SolverBase2D<cm_t, ExecutionBackend::CUDA> {
  using Base = SolverBase2D<cm_t, ExecutionBackend::CUDA>;

public:
  TMPCUDASolver2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {};

  ~TMPCUDASolver2D() = default;

  __host__ void solve(Lattice<2> &lattice,
                      const Params<2, cm_t> &params_) const override {
    double *fto, *ffrom, &norms;
    std::size_t allocation_size =
        lattice.grid.getArea() * sizeof(double) * D2Q9::ndir;

    // FIXME: add LBM_CUDA_CHECK
    LBM_CUDA_CHECK(cudaMalloc(&ffrom, allocation_size));
    LBM_CUDA_CHECK(cudaMalloc(&fto, allocation_size));
    LBM_CUDA_CHECK(cudaMalloc(&norms, allocation_size));

    // FIXME: Define the correct way to calculate blocks.

    // WARN: ExecutionContext is actually useless try to avoid it as much as
    // possible!!
    const dim3 blockDim = 1;
    const dim3 gridDim = 1;
    // FIXME: Needs to be instantiated.
    const cudaStream_t stream;

    // TODO: call init_equilibrium kernel

    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);

    for (auto iter = 0; iter < this->niters; iter++) {
      const bool save = iter % this->nskips == 0;
      update_stream_collide << gridDim, blockDim, 0,
          stream >> (lattice, cs, ffrom, fto, save);
      LBM_CUDA_CHECK(cudaGetLastError());
      std::swap(ffrom, fto);
      // TODO: Check whether all of this is really necessary
      if (save) {
        LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
        download_macroscopic(lattice, stream);
        write_norms(lattice);
      }
    }
  }

private:
  // TODO: check whether static is really needed
  static __global__ void
  init_equilibrium(const Lattice<2> &lattice,
                   std::vector<double> &part_stream) const {}

  // TODO: check whether static is really needed

  static __global__ void update_stream_collide(
      Lattice<2> &lattice, const CollisionStrategy<2, D2Q9, cm_t, CUDA> &cs,
      const double *__restrict__ ffrom, double *__restrict__ fto, int nx,
      int ny, double init_vel) {}

  static __device__ void field_index(int x, int y, std::size_t nx,
                                     std::size_t ndir) {
    // NOTE: More efficient memory access
    return ndir * (nx * y + x) + i;
  }

  static __device__ void apply_boundary_conditions(double *fto, double *ffrom,
                                                   int x, int y,
                                                   double localrho,
                                                   double init_vel) {
    // LEFT BOUNDARY: RESTING WALL
    if (x == 0) {
      fto[field_index(x, y, 1, 9)] = fto[field_index(x, y, 3, 9)];
      fto[field_index(x, y, 5, 9)] = fto[field_index(x, y, 7, 9)];
      fto[field_index(x, y, 8, 9)] = fto[field_index(x, y, 6, 9)];
    }

    // RIGHT BOUNDARY: RESTING WALL
    if (x == nx_ - 1) {
      fto[field_index(x, y, 3, 9)] = fto[field_index(x, y, 1, 9)];
      fto[field_index(x, y, 7, 9)] = fto[field_index(x, y, 5, 9)];
      fto[field_index(x, y, 6, 9)] = fto[field_index(x, y, 8, 9)];
    }

    // BOTTOM BOUNDARY: RESTING WALL
    if (y == 0) {
      fto[field_index(x, y, 4, 9)] = fto[field_index(x, y, 2, 9)];
      fto[field_index(x, y, 7, 9)] = fto[field_index(x, y, 5, 9)];
      fto[field_index(x, y, 8, 9)] = fto[field_index(x, y, 6, 9)];
    }

    // TOP BOUNDARY: MOVING WALL
    if (y == ny_ - 1) {
      // FIXME: Can't live with local encoding X*(((((
      fto[field_index(x, y, 2, 9)] =
          fto[field_index(x, y, 4, 9)] -
          2 * wi * localrho * (subx * ux + suby * uy) * 3;
      fto[field_index(x, y, 5, 9)] =
          fto[field_index(x, y, 7, 9)] -
          2 * wi * localrho * (subx * ux + suby * uy) * 3;
      fto[field_index(x, y, 6, 9)] =
          fto[field_index(x, y, 8, 9)] -
          2 * wi * localrho * (subx * ux + suby * uy) * 3;
    }
  }

  // FIXME: fix

  // Scarica la norma della velocità (già calcolata su device dal kernel
  // di collisione) e la inoltra agli observer registrati.
  void write_norms_from_device(cudaStream_t stream) const {
    std::vector<float> host_norms(area_);
    LBM_CUDA_CHECK(cudaMemcpyAsync(host_norms.data(), d_norm_,
                                   area_ * sizeof(float),
                                   cudaMemcpyDeviceToHost, stream));
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));

    std::vector<char> buf(host_norms.size() * sizeof(float));
    std::memcpy(buf.data(), host_norms.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }

  // Override richiesto da SolverBase2D (senza CUDASolver2D resta
  // astratta e non è istanziabile): i dati sono già sul device, la firma
  // con "lattice" serve solo a rispettare la virtuale della base, come fa
  // MPISolver2D::write_norms lato host tramite DataObservable.
  void write_norms(const Lattice<2> &lattice) const override {
    (void)lattice;
    write_norms_from_device(current_stream());
  }
};

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH
