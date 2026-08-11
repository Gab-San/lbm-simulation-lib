#ifndef __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH
#define __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH

// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

// CUDA
#include <cuda_runtime.h>

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {
namespace cuda_detail {

static __constant__ double c_wi[9];
static __constant__ int c_dirx[9];
static __constant__ int c_diry[9];
static __constant__ int c_opp[9];


// Costruisce i parametri device-side a partire dai Params host (una reference a
// Params<dim, cm_t>).
template <lbm::CollisionModel cm_t>
inline CudaCollisionParams<cm_t>
make_cuda_params(const lbm::Params<2, cm_t> &params) {
  if constexpr (cm_t == lbm::CollisionModel::BGK) {
    return CudaCollisionParams<cm_t>{params.tauinv, params.omtauinv};
  } else if constexpr (cm_t == lbm::CollisionModel::TRT) {
    return CudaCollisionParams<cm_t>{params.s_plus, params.s_minus};
  } else {
    static_assert(cm_t == lbm::CollisionModel::BGK,
                  "make_cuda_params: modello non implementato per CUDA.");
  }
}

// Indicizzazione identica a Grid<2>::field_index, ma su interi semplici
// (nessun oggetto Grid/std::vector è utilizzabile lato device ??? perchè??? da rivedere). 
// Layout -> "structure of arrays" (tutte le celle della direzione 0, poi tutte
// quelle della direzione 1, ...): per thread contigui che variano la sola
// x (come nei nostri kernel 2D) gli accessi restano coalescenti.
__device__ __forceinline__ std::size_t
d_field_index(int x, int y, int dir, int nx, int ny) {
  return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * dir +
         static_cast<std::size_t>(nx) * y + x;
}

// Collisione BGK, applicata da un singolo thread ad un singolo nodo
// (x, y). Il ciclo sulle 9 direzioni resta interno al thread: è un
// semplice for scritto per intero, nessuna direttiva di unrolling.
static __device__ __forceinline__ void
collide_bgk_node(double *__restrict__ f, int x, int y, int nx, int ny,
                 double rho, double ux, double uy,
                 const CudaCollisionParams<lbm::CollisionModel::BGK> &p) {
  const double u_sq = ux * ux + uy * uy;
  const double omusq = -1.5 * u_sq;

  const std::size_t plane =
      static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;

  for (int i = 0; i < 9; ++i) {
    const double cidotu = c_dirx[i] * ux + c_diry[i] * uy;
    const double feq =
        c_wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

    const std::size_t idx = plane * i + base;
    f[idx] = p.omtauinv * f[idx] + p.tauinv * feq;
  }
}


// Collisione TRT. Stessa logica dell'equivalente host: si processa ogni
// coppia (i, opp[i]) una sola volta (quando i <= opp[i]), la direzione di
// riposo (i == opp[i] == 0) è un caso a parte.
static __device__ __forceinline__ void
collide_trt_node(double *__restrict__ f, int x, int y, int nx, int ny,
                 double rho, double ux, double uy,
                 const CudaCollisionParams<lbm::CollisionModel::TRT> &p) {
  const double u_sq = ux * ux + uy * uy;
  const double c1 = -1.5 * u_sq;

  const std::size_t plane =
      static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  const std::size_t base = static_cast<std::size_t>(nx) * y + x;

  for (int i = 0; i < 9; ++i) {
    const int ib = c_opp[i];
    if (i > ib) {
      continue;
    }

    const double cidotu_i = c_dirx[i] * ux + c_diry[i] * uy;
    const double feq_i =
        c_wi[i] * rho * (1.0 + 3.0 * cidotu_i + 4.5 * cidotu_i * cidotu_i + c1);

    const std::size_t idx_i = plane * i + base;

    if (i == ib) {
      // Direzione di riposo: nessuna componente antisimmetrica.
      f[idx_i] = f[idx_i] - p.s_plus * (f[idx_i] - feq_i);
      continue;
    }

    const double cidotu_opp = c_dirx[ib] * ux + c_diry[ib] * uy;
    const double feq_opp =
        c_wi[ib] * rho *
        (1.0 + 3.0 * cidotu_opp + 4.5 * cidotu_opp * cidotu_opp + c1);

    const std::size_t idx_opp = plane * ib + base;

    const double fplus = 0.5 * (f[idx_i] + f[idx_opp]);
    const double fminus = 0.5 * (f[idx_i] - f[idx_opp]);
    const double fplus_eq = 0.5 * (feq_i + feq_opp);
    const double fminus_eq = 0.5 * (feq_i - feq_opp);

    f[idx_i] = f[idx_i] - p.s_plus * (fplus - fplus_eq) -
               p.s_minus * (fminus - fminus_eq);
    f[idx_opp] = f[idx_opp] - p.s_plus * (fplus - fplus_eq) +
                 p.s_minus * (fminus - fminus_eq);
  }
}


    // Dispatch a tempo di compilazione (CollisionStrategy::apply)
template <lbm::CollisionModel cm_t>
static __device__ __forceinline__ void
collide_node(double *__restrict__ f, int x, int y, int nx, int ny,
             double rho, double ux, double uy,
             const CudaCollisionParams<cm_t> &params) {
  if constexpr (cm_t == lbm::CollisionModel::BGK) {
    collide_bgk_node(f, x, y, nx, ny, rho, ux, uy, params);
  } else if constexpr (cm_t == lbm::CollisionModel::TRT) {
    collide_trt_node(f, x, y, nx, ny, rho, ux, uy, params);
  } else {
    static_assert(cm_t == lbm::CollisionModel::BGK,
                  "collide_node: modello non implementato per CUDA.");
  }
}

} //namespace cuda_detail
} //namespace lbm

#endif // __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH
