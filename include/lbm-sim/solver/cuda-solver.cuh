#ifndef __LBM_SIM_SOLVER_CUDA_SOLVER_CUH
#define __LBM_SIM_SOLVER_CUDA_SOLVER_CUH

#include "lbm-sim/backend.hpp"
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-operators/collision-strategy.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/constants.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/cuda/structs.cuh"
#include "lbm-sim/cuda/utils.cuh"
#include "lbm-sim/solver/solver-base.hpp"

// CUDA LIB
#include <cuda_runtime.h>

// C++ STD LIB
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
  const types::Coordinate<dim> p = cuda::thread_coordinate<dim>();

  // NOTE: Can this be replaced by an assertion?
  // Is this really needed?
  if (!grid.contains(p)) {
    return;
  }

  const double r = lattice_rho[grid.scalar_index(p)];
  const utils::Vector<double, dim> u = lattice_u[grid.scalar_index(p)];
  const double u_sq = utils::ops::dot(u, u);

  for (auto i = 0; i < VelocitySet::ndir; ++i) {
    const double cidotu = utils::ops::dot(cuda::vs_dir<dim, VelocitySet>[i], u);
    part_stream[grid.field_index(p, i, VelocitySet::ndir)] =
        cuda::vs_wi<VelocitySet>[i] * r *
        (1.0 + numbers::invcs_2 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
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

  const types::Coordinate<dim> p = cuda::thread_coordinate<dim>();

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
    const types::Coordinate<dim> src =
        p - cuda::vs_dir<dim, VelocitySet>[diridx];

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
  utils::Vector<double, dim> u;

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
    norms[s_idx] = __fsqrt_rn(static_cast<float>(utils::ops::dot(u, u)));
  }

  cs.apply(fp, p, r, u);

  // COPY BACK ON DEVICE BUFFER
  for (auto diridx = 0; diridx < VelocitySet::ndir; diridx++) {
    fto[grid.field_index(p, diridx, VelocitySet::ndir)] = fp[diridx];
  }
}

} // namespace cuda_detail

template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
class CUDASolver
    : public SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::CUDA> {
  using Base = SolverBase<dim, VelocitySet, cm_t, ExecutionBackend::CUDA>;

public:
  CUDASolver(const unsigned int iters_, const unsigned int frames_)
      : Base(iters_, frames_) {};

  ~CUDASolver() = default;

  __host__ void
  solve(Lattice<dim> &lattice,
        const CollisionParams<dim, cm_t> &params_) const override {

    cuda::upload_lattice_constants<dim, VelocitySet>();

    std::size_t area = lattice.grid.getArea();
    std::size_t allocation_size = area * sizeof(double) * VelocitySet::ndir;

    cuda::DeviceBuffer<double> ffrom(area * VelocitySet::ndir);
    cuda::DeviceBuffer<double> fto(area * VelocitySet::ndir);
    cuda::DeviceBuffer<types::boundary_t> d_boundary_mask(area);
    cuda::DeviceBuffer<float> norms(area);
    cuda::DeviceBuffer<double> d_rho(area);
    cuda::DeviceBuffer<utils::Vector<double, dim>> d_u(area);

    cuda::StreamHandler stream;

    // ----- DATA STRUCTURES INITIALIZATION -------
    d_rho.upload_async(lattice.rho, stream);
    d_u.upload_async(lattice.u, stream);
    d_boundary_mask.upload_async(lattice.boundary_mask, stream);
    // -----------------------------------------

    const dim3 block = cuda::create_block_of<dim>(dim == 2 ? 16 : 8);
    const dim3 grid_dims = cuda::ceil_div(lattice.grid.size, block);

    quill::Logger *solver_logger = logging::create_or_get_logger("solver");
    cuda::log_device_info(solver_logger);

    cuda_detail::init_equilibrium<dim, VelocitySet>
        <<<grid_dims, block, 0, stream>>>(ffrom.data(), d_rho.data(),
                                          d_u.data(), lattice.grid);

    LBM_CUDA_CHECK(cudaGetLastError());

    LOG_DEBUG(solver_logger, "Equilibrium Initialized...");

    for (auto iter = 0; iter < this->niters; iter++) {
      const bool save = iter % this->nskips == 0;
      const bool store_macroscopic = save || (iter + 1 == this->niters);

      cuda_detail::update_stream_collide<dim, VelocitySet, cm_t>
          <<<grid_dims, block, 0, stream>>>(
              ffrom.data(), fto.data(), d_boundary_mask.data(), norms.data(),
              d_rho.data(), d_u.data(), lattice.grid,
              CollisionStrategy<dim, VelocitySet, cm_t>(params_), lattice.pin,
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
  }

private:
  __host__ inline void download_macroscopic(
      const cuda::DeviceBuffer<double> &d_rho,
      const cuda::DeviceBuffer<utils::Vector<double, dim>> &d_u,
      Lattice<dim> &lattice, cudaStream_t stream) const {
    d_rho.download_async(lattice.rho, stream);
    d_u.download_async(lattice.u, stream);
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

  __host__ inline void write_norms(const cuda::DeviceBuffer<float> &d_norms,
                                   Grid<dim> grid, cudaStream_t stream) const {
    // Norm vector on host
    std::vector<float> h_norms(grid.getArea());
    d_norms.download_async(h_norms, stream);
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));

    std::vector<char> buf(h_norms.size() * sizeof(float));
    std::memcpy(buf.data(), h_norms.data(), buf.size());
    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm

#endif // __LBM_SIM_SOLVER_CUDA_SOLVER_CUH
