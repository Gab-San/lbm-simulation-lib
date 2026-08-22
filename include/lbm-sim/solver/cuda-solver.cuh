#ifndef __LBM_SIM_SOLVER_SOLVER_CUDA_CUH
#define __LBM_SIM_SOLVER_SOLVER_CUDA_CUH

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/boundaries.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/backend.hpp"
#include "lbm-sim/cuda/utils.cuh"

// CUDA LIB
#include <cuda_runtime.h>

// C++ STD LIB
#include <cmath>
#include <cstring>
#include <utility>
#include <vector>

namespace lbm {
namespace cuda_detail {

template <types::dim_t dim, typename VelocitySet>
__global__ void
init_equilibrium(double *__restrict__ part_stream,
                 const double *__restrict__ lattice_rho,
                 const utils::Vector<double, dim> *__restrict__ lattice_u,
                 const Grid<dim> grid) {
  const types::Coordinate<dim> p(blockIdx.x * blockDim.x + threadIdx.x,
                                 blockIdx.y * blockDim.y + threadIdx.y);

  if (!grid.contains(p)) {
    return;
  }

  const double r = lattice_rho[grid.scalar_index(p)];
  const utils::Vector<double, 2> u = lattice_u[grid.scalar_index(p)];
  const double u_sq = utils::ops::dot(u, u);

  for (auto i = 0; i < VelocitySet::ndir; ++i) {
    const double cidotu = utils::ops::dot(cuda::vs_dir<dim, VelocitySet>[i], u);
    part_stream[grid.field_index(p, i, VelocitySet::ndir)] =
        cuda::vs_wi<VelocitySet>[i] * r *
        (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
  }
}

template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
__global__ void update_stream_collide(
    const double *__restrict__ const ffrom, double *__restrict__ fto,
    const types::boundary_t *__restrict__ boundary_mask,
    float *__restrict__ norms, double *__restrict__ lattice_rho,
    utils::Vector<double, dim> *__restrict__ lattice_u, const Grid<dim> grid,
    const CollisionStrategy<dim, VelocitySet, cm_t> cs, const double pin,
    const double pout, const bool store_macroscopic) {
  double fp[VelocitySet::ndir];

  // NOTE: Can this be dimension agnostic?
  const types::Coordinate<2> p(blockIdx.x * blockDim.x + threadIdx.x,
                               blockIdx.y * blockDim.y + threadIdx.y);

  // NOTE: Can this be replaced by an assertion?
  // Is this really needed?
  if (!grid.contains(p)) {
    return;
  }

  double r_wall = 0.0;
  for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
    r_wall += ffrom[grid.field_index(p, diridx, VelocitySet::ndir)];
  }

  // STREAMING + HALFWAY COLLISION
  for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
    const types::Coordinate<2> src = p - cuda::vs_dir<dim, VelocitySet>[diridx];

    if (!grid.contains(src)) {
      // if source node is external it is on a boundary node
      Solid::apply_boundary_condition<dim, VelocitySet>(
          fp, ffrom, diridx, grid, boundary_mask, lattice_rho, lattice_u, p,
          r_wall, cs.params.init_vel, pin, pout);
    } else {
      // if source node is internal stream it.
      fp[diridx] = ffrom[grid.field_index(src, diridx, VelocitySet::ndir)];
    }
  }

  // COMPUTE MACROSCOPIC VARIABLES
  // rho = sum_i fi
  // rho*u = sum_i fi * ci

  double r = 0.0;
  utils::Vector<double, 2> u(0, 0);

  for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
    r += fp[diridx];
    u += cuda::vs_dir<dim, VelocitySet>[diridx] * fp[diridx];
  }

  // u = (sum_i fi * ci) / rho
  u /= r;

  // STORE computed macroscopic values

  // Store macroscopic fields for requested frames and for the final state.
  const auto s_idx = grid.scalar_index(p);
  if (store_macroscopic ||
      boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_INLET ||
      boundary_mask[s_idx] == Solid::PRESSURE_PERIODIC_OUTLET) {
    lattice_rho[s_idx] = r;
    lattice_u[s_idx] = u;
  }

  if (store_macroscopic) {
    norms[s_idx] = static_cast<float>(std::sqrt(utils::ops::dot(u, u)));
  }

  cs.apply(fp, p, r, u);

  // COPY BACK ON DEVICE BUFFER
  for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
    fto[grid.field_index(p, diridx, VelocitySet::ndir)] = fp[diridx];
  }
}

} // namespace cuda_detail

template <enum CollisionModel cm_t>
class CUDASolver2D : public SolverBase2D<cm_t, ExecutionBackend::CUDA> {
  using Base = SolverBase2D<cm_t, ExecutionBackend::CUDA>;

public:
  CUDASolver2D(const unsigned int num_iters_, const unsigned int num_frames_)
      : Base(num_iters_, num_frames_) {};

  ~CUDASolver2D() = default;

  __host__ void solve(Lattice<2> &lattice,
                      const CollisionParams<2, cm_t> &params_) const override {

    cuda::upload_lattice_constants<2, D2Q9>();
    double *fto, *ffrom, *d_rho;
    types::boundary_t *d_boundary_mask;
    float *norms;
    utils::Vector<double, 2> *d_u;

    cudaStream_t stream;
    LBM_CUDA_CHECK(cudaStreamCreate(&stream));

    std::size_t area = lattice.grid.getArea();
    std::size_t allocation_size = area * sizeof(double) * D2Q9::ndir;

    LBM_CUDA_CHECK(cudaMalloc(&ffrom, allocation_size));
    LBM_CUDA_CHECK(cudaMalloc(&fto, allocation_size));
    LBM_CUDA_CHECK(
        cudaMalloc(&d_boundary_mask, area * sizeof(types::boundary_t)));

    LBM_CUDA_CHECK(cudaMalloc(&norms, area * sizeof(float)));
    LBM_CUDA_CHECK(cudaMalloc(&d_rho, area * sizeof(double)));
    LBM_CUDA_CHECK(cudaMalloc(&d_u, area * sizeof(utils::Vector<double, 2>)));

    LBM_CUDA_CHECK(cudaMemcpyAsync(d_rho, lattice.rho.data(),
                                   area * sizeof(double),
                                   cudaMemcpyHostToDevice, stream));
    LBM_CUDA_CHECK(cudaMemcpyAsync(d_u, lattice.u.data(),
                                   area * sizeof(utils::Vector<double, 2>),
                                   cudaMemcpyHostToDevice, stream));
    LBM_CUDA_CHECK(cudaMemcpyAsync(
        d_boundary_mask, lattice.boundary_mask.data(),
        area * sizeof(types::boundary_t), cudaMemcpyHostToDevice, stream));

    // 16x16 is a conservative 2D default. ceil_div below guarantees complete
    // domain coverage; device-specific tuning can be added independently.
    const dim3 blockDim(16, 16);
    const dim3 gridDim(cuda_detail::ceil_div(lattice.grid.size.x, blockDim.x),
                       cuda_detail::ceil_div(lattice.grid.size.y, blockDim.y));

    cuda_detail::init_equilibrium<2, D2Q9>
        <<<gridDim, blockDim, 0, stream>>>(ffrom, d_rho, d_u, lattice.grid);

    for (auto iter = 0; iter < this->niters; iter++) {
      const bool save = iter % this->nskips == 0;
      const bool store_macroscopic = save || (iter + 1 == this->niters);

      cuda_detail::update_stream_collide<2, D2Q9, cm_t>
          <<<gridDim, blockDim, 0, stream>>>(
              ffrom, fto, d_boundary_mask, norms, d_rho, d_u, lattice.grid,
              CollisionStrategy<2, D2Q9, cm_t>(params_), lattice.pin,
              lattice.pout, store_macroscopic);

      LBM_CUDA_CHECK(cudaGetLastError());
      std::swap(ffrom, fto);

      if (save) {
        LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
        download_macroscopic(d_rho, d_u, lattice, stream);
        write_norms(norms, lattice.grid, stream);
      }
    }

    // The final iteration always stores rho/u even when it is not a frame.
    download_macroscopic(d_rho, d_u, lattice, stream);
    LBM_CUDA_CHECK(cudaStreamDestroy(stream));

    LBM_CUDA_CHECK(cudaFree(ffrom));
    LBM_CUDA_CHECK(cudaFree(fto));
    LBM_CUDA_CHECK(cudaFree(d_boundary_mask));
    LBM_CUDA_CHECK(cudaFree(norms));
    LBM_CUDA_CHECK(cudaFree(d_rho));
    LBM_CUDA_CHECK(cudaFree(d_u));
  }

private:
  __host__ void download_macroscopic(const double *const d_rho,
                                     const utils::Vector<double, 2> *const d_u,
                                     Lattice<2> &lattice,
                                     cudaStream_t stream) const {
    const auto area = lattice.grid.getArea();

    LBM_CUDA_CHECK(cudaMemcpyAsync(lattice.rho.data(), d_rho,
                                   area * sizeof(double),
                                   cudaMemcpyDeviceToHost, stream));

    LBM_CUDA_CHECK(cudaMemcpyAsync(lattice.u.data(), d_u,
                                   area * sizeof(utils::Vector<double, 2>),
                                   cudaMemcpyDeviceToHost, stream));

    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  __host__ inline void write_norms(float *dnorms, Grid<2> grid,
                                   cudaStream_t stream) const {
    // Norm vector on host
    std::vector<float> h_norms(grid.getArea());
    LBM_CUDA_CHECK(cudaMemcpyAsync(h_norms.data(), dnorms,
                                   h_norms.size() * sizeof(float),
                                   cudaMemcpyDeviceToHost, stream));
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));

    std::vector<char> buf(h_norms.size() * sizeof(float));
    std::memcpy(buf.data(), h_norms.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH
