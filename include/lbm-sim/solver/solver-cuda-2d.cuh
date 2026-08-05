#ifndef __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH
#define __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH


// NOTA IMPORTANTE SULLA BUILD:
// va incluso da un file .cu o compilato con nvcc NON con g++/clang


// LBM SIM LIB
#include "lbm-sim/solver/solver-base.hpp"
#include "lbm-sim/solver/solver-2d.hpp" // needed AsyncBinaryWriter
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/collision-operators/collision-strategy-cuda.cuh"
#include "lbm-sim/backend/metadata.hpp"

// COLLISION DETECTION LIB
#include "collision-detection/core/types.hpp"
#include "collision-detection/core/operators.hpp" 

// CUDA
#include <cuda_runtime.h>

// C++ STANDARD LIB
#include <cmath>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include <assert.h>

namespace lbm {

// Error checks for async errors
#define LBM_CUDA_CHECK(expr)
  do {
    const cudaError_t _lbm_err = (expr);
    if (_lbm_err != cudaSuccess) {
      _lbm_oss << "CUDA ERROR " <<  cudaGetErrorString(err) << " at " 
      << __FILE__ << " : " << __LINE__;
      throw std::runtime_error(_lbm_oss.str());
    }
  }
  while (0)

namespace cuda_detail{

  inline unsigned int ceil_div(unsigned int a, unsigned int b){
    return ( a + b - 1 ) / b;
  }

  struct BoundaryNode {
  int x;
  int y;
  int is_moving;
  };
  /* boundary node serve per gestire indice appiattito, sostituisce sul device la struttura obstacles[]
  /get perimeter usata da host che dipende da collision-detection library. (da vedere se funzia???)
  costruita e caricata una sola volta su device, in costruzione del solver, la geom di ostacoli e pareti
  non cambia tra iter e l'altra, cambia solo la velox della parete mobile */


// KERNEL 1: streaming + calc densità temporanea
// uso un thread per nodo, lavora su dati dello stesso "worker" (ha appena scritto le 9 direz
// del proprio nodo) NO DIPENDENZA CROSS-THREAD tra kernel stream e kernel density
// lancio un kernel in meno. 

// caso sorgente fuori dal dominio -> se p - dir[i] cade fuori da griglia, lo slot fto NON viene
// scritto (si legge il val già presente, in openMp era continue). ogni nodo il cui vicino esce 
// dal dominio sia registrato come nodo di bordo/ostacolo + corretto dal kernel boundary condition

static __global__ void kernel_stream_density(const double *__restrict__ ffrom,
                                        double *__restrict__ fto, int nx, int ny) {
  const int x = blockIdx.x * blockDim.x + threadIdx.x; //offset cuda
  const int y = blockIdx.y * blockDim.y + threadIdx.y;
  if (x >= nx || y >= ny){
      // Magari qui si potrebbe lanciare un errore.
      // Sarebbe da vedere, perchè non so se si possa limitare la grandezza
      // del kernel in base alla griglia (non ricordo se sia dinamico o
      // statico).
    return;
  }

    // Qui si potrebbe usare il metodo getArea() di grid (forse non è il massimo
    // passare l'intera grid ci dovrei pensare)
    // Oppure possiamo definire un metodo per poter calcolare il prodotto tra
    // componenti di un vettore e calcolare così la grandezza di tutto.
  const std::size_t plane = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;

  double r = 0.0;

    // fto e fdst qui dovrebbero essere o grandi
    // quanto l'intera porzione di griglia di cui il blocco CUDA
    // si occupa oppure si potrebbero usare dei semplici array grossi
    // ndir
  for (int i = 0; i < 9; ++i) {
    const int sx = x - c_dirx[i];
    const int sy = y - c_diry[i];
    // Sarebbe meglio usare il calcolo vettoriale.
    const std::size_t dst = plane * i + base;
    if (sx >= 0 && sx < nx && sy >= 0 && sy < ny) {
      const double fi =
          ffrom[plane * i + static_cast<std::size_t>(nx) * sy + sx];
      fto[dst] = fi;
      r += fi;
    } else {
      // continue, se serve fa update in boundary conditions
      r += fto[dst];
    }
  }

  // Dipende in questo caso rho_tmp a cosa ti serve/
    // Nel codice nuovo rho_tmp serve solamente per applicare le condizioni al
    // contorno. L'applicazione di queste verrà probabilmente spostata durante
    // lo streaming. Questo rende meno problematico il caso delle condizioni al
    // contorno.

  rho_tmp[base] = r;
}

  // PER I VARI AGGIORNAMENTI AL CODICE (che ne semplificano la scrittura)
    // guardare
    // il branch optimizations/...


// KERNEL 2: Boundary conditions (bounce back semplice / con correzione per parete mobile)
// Un thread per nodo di bordo, non per nodo di griglia: la lista è più corta del dominio 
// completo. usa grid 1D dedicato senza scartare al kernel domain wide il lavoro sui nodi 
// che non sono di bordo


static __global__ void boundary_kernel(double *__restrict__ f, const double *__restrict__ rho_tmp,
                const BoundaryNode *__restrict__ boundary, int n_boundary,
                int nx, int ny, double ux0, double uy0) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n_boundary){
    return;
  }

  const BoundaryNode node = boundary[idx];
  const int x = node.x;
  const int y = node.y;

  #if !defined(NDEBUG)
    assert(x >= 0 && x < nx && y >= 0 && y < ny);
  #endif

  const std::size_t plane = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;
  const double r = rho_tmp[base];

  // La direzione 0 (riposo) viene saltata: coincide col nodo ostacolo
  // stesso, esattamente come nella versione host.
  for (int i = 1; i < 9; ++i) {
    const int tx = x + c_dirx[i];
    const int ty = y + c_diry[i];
    if (tx < 0 || tx >= nx || ty < 0 || ty >= ny) {
      continue;
    }

    const int ib = c_opp[i];
    const std::size_t f_i = plane * i + base;
    const std::size_t f_ib = plane * ib + base;

    if (node.is_moving) {
      const double cidotu = c_dirx[i] * ux0 + c_diry[i] * uy0;
      // c_s = 1/sqrt(3) => 1/c_s^2 = 3:  FIXME di prima
      // parametro invece di una costante
      f[f_i] = f[f_ib] - 2.0 * c_wi[i] * r * cidotu * 3.0;
    } else {
      f[f_i] = f[f_ib];
    }
  }
}


// KERNEL 3: calcolo delle quantità macroscopiche + collisione, lo stesso 
// thread che ha appena sommato rho/u esegue subito collisione, riusando ux/uy
// dai registri senza ripassare da global memory!!
// Quando save == true, oltre a rho/u viene scritta anche la norma della
// velocità in norm_out (write norms)
// Così l'intera pipeline di calcolo resta esclusivamente su GPU e la CPU
// si limita a spostare su disco pochi float già pronti.

template <lbm::CollisionModel cm_t>
static __global__ void kernel_macroscopic_and_collide(
    double *__restrict__ f, double *__restrict__ rho_out,
    double2 *__restrict__ u_out, float *__restrict__ norm_out, int nx,
    int ny, bool save, CudaCollisionParams<cm_t> params) {
  const int x = blockIdx.x * blockDim.x + threadIdx.x;
  const int y = blockIdx.y * blockDim.y + threadIdx.y;
  if (x >= nx || y >= ny) {
    return;
  }

  const std::size_t plane = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;

  double rho = 0.0;
  double ux = 0.0;
  double uy = 0.0;

  for (int i = 0; i < 9; ++i) {
    const double fi = f[plane * i + base];
    rho += fi;
    ux += c_dirx[i] * fi;
    uy += c_diry[i] * fi;
  }
  ux /= rho;
  uy /= rho;

  if (save) {
    rho_out[base] = rho;
    u_out[base] = make_double2(ux, uy);
    norm_out[base] = static_cast<float>(sqrt(ux * ux + uy * uy));
  }

  collide_node<cm_t>(f, x, y, nx, ny, rho, ux, uy, params); //def in collision-strategy-cuda
}


// KERNEL 0: inizializzazione all'equilibrio

static __global__ void
kernel_init_equilibrium(double *__restrict__ f, const double *__restrict__ rho,
                        const double2 *__restrict__ u, int nx, int ny) {
  const int x = blockIdx.x * blockDim.x + threadIdx.x;
  const int y = blockIdx.y * blockDim.y + threadIdx.y;
  if (x >= nx || y >= ny) {
    return;
  }

  const std::size_t plane = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;

  const double r = rho[base];
  const double2 vel = u[base];
  const double u_sq = vel.x * vel.x + vel.y * vel.y;

  for (int i = 0; i < 9; ++i) {
    const double cidotu = c_dirx[i] * vel.x + c_diry[i] * vel.y;
    f[plane * i + base] = c_wi[i] * r * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * u_sq);
  }
}

// Copia le costanti del set di velocità D2Q9 (definite host-side in velocity-sets.hpp) 
// nella constant memory dichiarata in collision-strategy-cuda.cuh. Va eseguita una volta 
// per processo, prima di lanciare qualunque kernel: viene chiamata dal costruttore del
// solver.
inline void upload_lattice_constants() {
  double h_wi[9];
  int h_dirx[9];
  int h_diry[9];
  int h_opp[9];

  for (int i = 0; i < 9; ++i) {
    h_wi[i] = D2Q9::wi[i];
    h_dirx[i] = D2Q9::dir[i].dx;
    h_diry[i] = D2Q9::dir[i].dy;
    h_opp[i] = D2Q9::opp[i];
  }

  LBM_CUDA_CHECK(cudaMemcpyToSymbol(c_wi, h_wi, sizeof(h_wi)));
  LBM_CUDA_CHECK(cudaMemcpyToSymbol(c_dirx, h_dirx, sizeof(h_dirx)));
  LBM_CUDA_CHECK(cudaMemcpyToSymbol(c_diry, h_diry, sizeof(h_diry)));
  LBM_CUDA_CHECK(cudaMemcpyToSymbol(c_opp, h_opp, sizeof(h_opp)));
}

}//namespace cuda_detail



// CUDASolver2D:
//parallelizzazione: un thread cuda per nodo di griglia (x,y), mappato su griglia 2D (FIX ME GENERALIZZARE)
//di blocchi/thread. il ciclo sulle 9 direz resta sequenziale al singolo thread. in collisione
//e non c'è lavoro condiviso tra thread diversi per direzioni e tenerle nello stesso thread permette
//di accumulare rho/ux/uy in registri invece di global mem. 
//BC usano un grid 1D dedicato ai soli thread non interni.
//f struct of array,  (stride nx*ny tra direz), identico a Grid<2>::field_index: per un warp che varia la
// sola x, questo produce coalesced accesses. (FIXME DA VERIFICARE)
template <enum CollisionModel cm_t>
class CUDASolver2D : public SolverBase2D<cm_t, ExecutionBackend::CUDA> {
  using Base = SolverBase2D<cm_t, ExecutionBackend::CUDA>;

public:
  CUDASolver2D(const unsigned int num_iters_, const unsigned int num_frames_,
              const Structure<2> &strt_, const std::string &out_path,
              const ExecutionContext<ExecutionBackend::CUDA> &ctx = {})
      : Base(num_iters_, num_frames_, strt_), norms_writer(out_path),
        ctx_(ctx) {
    cuda_detail::upload_lattice_constants();
    upload_boundary_nodes();
  }

  ~CUDASolver2D() override {
    free_field_buffers();
    if (d_boundary_ != nullptr) {
      cudaFree(d_boundary_);
    }
  }

  CUDASolver2D(const CUDASolver2D &) = delete;
  CUDASolver2D &operator=(const CUDASolver2D &) = delete;

  void init_equilibrium(Grid<2> &grid,
                        std::vector<double> &part_stream) const override {
    ensure_device_buffers(grid);
    const cudaStream_t stream = current_stream();

    LBM_CUDA_CHECK(cudaMemcpyAsync(d_rho_, grid.rho.data(),
                                   area_ * sizeof(double),
                                   cudaMemcpyHostToDevice, stream));

    std::vector<double2> h_u(area_);
    for (std::size_t k = 0; k < area_; ++k) {
      h_u[k] = make_double2(grid.u[k].dx, grid.u[k].dy);
    }
    LBM_CUDA_CHECK(cudaMemcpyAsync(d_u_, h_u.data(), area_ * sizeof(double2),
                                   cudaMemcpyHostToDevice, stream));

    const dim3 block(ctx_.block_x, ctx_.block_y);
    const dim3 gridDim(cuda_detail::ceil_div(static_cast<unsigned int>(nx_),
                                             block.x),
                       cuda_detail::ceil_div(static_cast<unsigned int>(ny_),
                                             block.y));

    cuda_detail::kernel_init_equilibrium<<<gridDim, block, 0, stream>>>(
        d_f_a_, d_rho_, d_u_, nx_, ny_);
    LBM_CUDA_CHECK(cudaGetLastError());

    part_stream.resize(area_ * D2Q9::ndir);
    LBM_CUDA_CHECK(cudaMemcpyAsync(part_stream.data(), d_f_a_,
                                   area_ * D2Q9::ndir * sizeof(double),
                                   cudaMemcpyDeviceToHost, stream));
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
  }

   void solve(Grid<2> &grid, const Params<2, cm_t> &params_,
             std::vector<double> &ffrom,
             std::vector<double> &fto) const override {
    ensure_device_buffers(grid);
    const cudaStream_t stream = current_stream();
    const std::size_t fsize = area_ * D2Q9::ndir;

    // ffrom fto come in host
    LBM_CUDA_CHECK(cudaMemcpyAsync(d_f_a_, ffrom.data(), fsize * sizeof(double),
                                   cudaMemcpyHostToDevice, stream));
    LBM_CUDA_CHECK(cudaMemcpyAsync(d_f_b_, fto.data(), fsize * sizeof(double),
                                   cudaMemcpyHostToDevice, stream));

    const cuda_detail::CudaCollisionParams<cm_t> cparams =
        cuda_detail::make_cuda_params(params_);
    const double ux0 = params_.init_vel.dx;
    const double uy0 = params_.init_vel.dy;

    write_header(grid);

    double *d_from = d_f_a_;
    double *d_to = d_f_b_;

    const dim3 block(ctx_.block_x, ctx_.block_y);
    const dim3 gridDim(cuda_detail::ceil_div(static_cast<unsigned int>(nx_),
                                             block.x),
                       cuda_detail::ceil_div(static_cast<unsigned int>(ny_),
                                             block.y));

    const unsigned int bc_threads = 128;
    const dim3 blockB(bc_threads);
    const dim3 gridB(cuda_detail::ceil_div(
        static_cast<unsigned int>(std::max(n_boundary_, 1)), bc_threads));

    for (unsigned int iter = 0; iter < this->niters; ++iter) {
      const bool save = (this->nskips > 0) && (iter % this->nskips == 0);

      cuda_detail::kernel_stream_and_density<<<gridDim, block, 0, stream>>>(
          d_from, d_to, d_rho_tmp_, nx_, ny_);

      if (n_boundary_ > 0) {
        cuda_detail::kernel_boundary_conditions<<<gridB, blockB, 0, stream>>>(
            d_to, d_rho_tmp_, d_boundary_, n_boundary_, nx_, ny_, ux0, uy0);
      }

      cuda_detail::kernel_macroscopic_and_collide<cm_t>
          <<<gridDim, block, 0, stream>>>(d_to, save ? d_rho_ : nullptr,
                                          save ? d_u_ : nullptr,
                                          save ? d_norm_ : nullptr, nx_, ny_,
                                          save, cparams);
      LBM_CUDA_CHECK(cudaGetLastError());

      std::swap(d_from, d_to);

      if (save) {
        LBM_CUDA_CHECK(cudaStreamSynchronize(stream));
        download_macroscopic(grid, stream);
        write_norms(stream);
      }
    }

    LBM_CUDA_CHECK(cudaMemcpyAsync(ffrom.data(), d_from, fsize * sizeof(double),
                                   cudaMemcpyDeviceToHost, stream));
    LBM_CUDA_CHECK(cudaMemcpyAsync(fto.data(), d_to, fsize * sizeof(double),
                                   cudaMemcpyDeviceToHost, stream));
    LBM_CUDA_CHECK(cudaStreamSynchronize(stream));

private:
  //nodi di bordo
  //current stream
  //ensure device buffers 
  //free buffers
  //upload boundary nodes
  //save macroscopic quantities
  //write norms
  

}

} //namespace lbm

#endif // __LBM_SIM_SOLVER_SOLVER_2D_CUDA_CUH