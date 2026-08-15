#ifndef __LBM_SIM_SOLVER_SOLVER_CUDA_CUH
#define __LBM_SIM_SOLVER_SOLVER_CUDA_CUH

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

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
    double *fto, *ffrom, *norms;
    std::size_t allocation_size =
        lattice.grid.getArea() * sizeof(double) * D2Q9::ndir;

    // FIXME: add LBM_CUDA_CHECK
    LBM_CUDA_CHECK(cudaMalloc(&ffrom, allocation_size));
    LBM_CUDA_CHECK(cudaMalloc(&fto, allocation_size));
    LBM_CUDA_CHECK(cudaMalloc(&norms, allocation_size));

    // FIXME: Define the correct way to calculate blocks.

    // WARN: ExecutionContext is actually useless try to avoid it as much as
    // possible!!

    // NOTE: Is there any way to make this dimension independent? (Maybe with an
    // helper function)
    const dim3 blockDim(16, 16);
    const dim3 gridDim(cuda_detail::ceil_div(lattice.grid.size.x, blockDim.x), cuda_detail::ceil_div(lattice.grid.size.y, blockDim.y);

    // FIXME: Needs to be instantiated.
    const cudaStream_t stream;

    // TODO: call init_equilibrium kernel

    const CollisionStrategy<2, D2Q9, cm_t, OPEN_MP> cs(params_);

    for (auto iter = 0; iter < this->niters; iter++) {
      const bool save = iter % this->nskips == 0;
      // FIXME: Add stream (if needed SPOILER: PROBABILMENTE SI)
      update_stream_collide<<<gridDim, blockDim>>>(lattice, cs, ffrom, fto,
                                                   norms, save);
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

  // TODO: check whether lattice is viably passable (otherwise it might be
  // adapted)
  static __global__ void update_stream_collide(
      Lattice<2> &lattice, const CollisionStrategy<2, D2Q9, cm_t, CUDA> &cs,
      const double *__restrict__ ffrom, double *__restrict__ fto,
      double *__restrict__ norms, double init_vel) {
    __shared__ double *bfto, bffrom;
    const types::Coordinate<2> p(blockIdx.x + blockDim.x + threadIdx.x,
                                 blockIdx.y * blockDim.y + threadIdx.y);

    // NOTE: Can this be replaced by an assertion?
    // Is this really needed?
    if (!lattice.grid.contains(p)) {
      return;
    }

    // FIXME: Access to shared memory is completely wrong!
    // Shared memory is never allocated!!

    // Copy on block local memory
    // NOTE: Is this really necessary? Does it speed up computation?
    double r_wall = 0.0;
    for (auto diridx = 0; diridx < D2Q9::ndir; diridx++) {
      bfto[lattice.grid.field_index(p, diridx, D2Q9::ndir)] =
          fto[lattice.grid.field_index(p, diridx, D2Q9::ndir)];
      bffrom[lattice.grid.field_index(p, dir, D2Q9::ndir)] =
          ffrom[lattice.grid.field_index(p, dir, D2Q9::ndir)];
      r_wall += bffrom[lattice.grid.field_index(p, dir, D2Q9::ndir)];
    }

    // STREAMING + HALFWAY COLLISION
    for (auto diridx = 0; diridx < D2Q9::ndir; diridx++) {
      const types::Coordinate<2> src = p - D2Q9::dir[diridx];
      if (!lattice.grid.contains(src)) {
        // if source node is external it is on a boundary node
        apply_boundary_conditions(bfto, bffrom, p.x, p.y, r_wall,
                                  cs.params.init_vel);
      } else {
        // if source node is internal stream it.
        bfto[lattice.grid.field_index(p, dir, D2Q9::ndir)] =
            bffrom[lattice.grid.field_index(src, diridx, D2Q9::ndir)];
      }
    }

    // COMPUTE MACROSCOPIC VARIABLES
    // rho = sum_i fi
    // rho*u = sum_i fi * ci

    double r = 0.0;
    utils::Vector<double, 2> u(0, 0);

    for (auto diridx = 0; diridx < D2Q9::ndir; diridx++) {
      const double fi = bfto[lattice.grid.field_index(p, diridx, D2Q9::ndir)];
      rho += fi;
      u += D2Q9::dir[i] * fi;
    }

    // u = (sum_i fi * ci) / rho
    u /= r;

    // STORE computed macroscopic values

    // FIXME: STORING NEED TO BE CORRECTED
    // Lattice probably is not compatible with cuda (can it be made compatible?)
    if (save) {
      const unsigned int s_idx = lattice.grid.scalar_index(p);
      // NOTE: use device buffers and then copy back on lattice?
      lattice.rho[s_idx] = r;
      lattice.u[s_idx] = u;
      norms = static_cast<float>(std::sqrt(utils::ops::dot(u, u)));
    }

    // FIXME: This will not work
    collide_node<cm_t>(f, x, y, nx, ny, rho, ux, uy,
                       params); // def in collision-strategy-cuda

    // COPY BACK ON DEVICE BUFFER
    for (auto diridx = 0; diridx < D2Q9::ndir; diridx++) {
      fto[lattice.grid.field_index(p, diridx, D2Q9::ndir)] =
          bfto[lattice.grid.field_index(p, diridx, D2Q9::ndir)];
    }
  }

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

  // FIXME: This functions do not work.
  // They must be adapted or the code must be corrected in order to use these
  // functions

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
  void write_norms(const Lattice<2> &lattice) const {
    (void)lattice;
    write_norms_from_device(current_stream());
  }
};

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH
