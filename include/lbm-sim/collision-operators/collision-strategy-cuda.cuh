#ifndef __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH
#define __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH

// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

// CUDA
#include <cuda_runtime.h>

// C++ STANDARD LIB
#include <cstddef>

namespace lbm {
namespace cuda_detail {
//
// //Dichiarazione di direzioni/pesi/opposti della griglia (visibile a tutti i
// thread device,
// //con cache dedicata). Popolati una sola volta, lato host, da
// //upload_lattice_constants() (chiamata da costruttore da CUDASolver2D).
// __constant__ int c_dirx[9];
// __constant__ int c_diry[9];
// __constant__ double c_wi[9];
// __constant__ int c_opp[9]; // indice della direzione opposta (per il
// bounce-back)
//
// //carica le costanti della griglia in __constant__ memory a partire dai
// valori host di
// //D2Q9::dir/wi. l'opposto di ogni direzione è calcolato per confronto diretto
// sui vettori
// //-> non serve conoscere/duplicare la convenzione di ordinamento usata da
// D2Q9. inline void upload_lattice_constants() {
//   static bool uploaded = false;
//   if (uploaded) {
//     return;
//   }
//
//   int h_dirx[D2Q9::ndir];
//   int h_diry[D2Q9::ndir];
//   double h_wi[D2Q9::ndir];
//   int h_opp[D2Q9::ndir];
//
//   for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
//     h_dirx[i] = static_cast<int>(D2Q9::dir[i].dx);
//     h_diry[i] = static_cast<int>(D2Q9::dir[i].dy);
//     h_wi[i] = D2Q9::wi[i];
//   }
//
//   for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
//     int opp = 0;
//     for (unsigned int j = 0; j < D2Q9::ndir; ++j) {
//       if (h_dirx[j] == -h_dirx[i] && h_diry[j] == -h_diry[i]) {
//         opp = static_cast<int>(j);
//         break;
//       }
//     }
//     h_opp[i] = opp;
//   }
//
//   cudaMemcpyToSymbol(c_dirx, h_dirx, sizeof(h_dirx));
//   cudaMemcpyToSymbol(c_diry, h_diry, sizeof(h_diry));
//   cudaMemcpyToSymbol(c_wi, h_wi, sizeof(h_wi));
//   cudaMemcpyToSymbol(c_opp, h_opp, sizeof(h_opp));
//
//   uploaded = true;
// }
}

namespace cuda {

template <types::dim_t dim, typename VelocitySet>
static __device__ __forceinline__ void
collide_bgk_node(double *__restrict__ fp, const types::Coordinate<dim> p,
                 const double rho, const utils::Vector<double, dim> u,
                 const Params<dim, CollisionModel::BGK> &params) {
  using utils::ops::dot;
  const double omusq = -1.5 * dot(u, u);

  for (int diridx = 0; diridx < VelocitySet::ndir; ++diridx) {
    const double cidotu = dot(vs_dir<dim, VelocitySet>[diridx], u);
    const double feq = vs_wi<VelocitySet>[diridx] * rho *
                       (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

    // RELAX TO EQUILIBRIUM
    fp[diridx] = params.omtauinv * fp[diridx] + params.tauinv * feq;
  }
}

template <types::dim_t dim, typename VelocitySet>
static __device__ __forceinline__ void
collide_trt_node(double *__restrict__ f, const types::Coordinate<dim> p,
                 const double rho, const utils::Vector<double, dim> u,
                 const Params<dim, CollisionModel::TRT> &params) {
  using utils::ops::dot;

  const double omusq = -1.5 * dot(u, u);

  for (auto i = 0; i < VelocitySet::ndir; ++i) {
    const int iopp = vs_opp<VelocitySet>[i];

    if (i > iopp) {
      continue;
    }

    const double cidotu_i = dot(vs_dir<dim, VelocitySet>[i], u);
    const double feq_i =
        vs_wi<VelocitySet>[i] * rho *
        (1.0 + 3.0 * cidotu_i + 4.5 * cidotu_i * cidotu_i + omusq);

    if (i == iopp) {
      // TODO: THIS COMMENT IS SHIT ENGLISH
      //
      // Center Direction: no antisymmetric component.
      f[i] = f[i] - params.s_plus * (f[i] - feq_i);
      continue;
    }

    const double cidotu_opp = dot(vs_dir<dim, VelocitySet>[iopp], u);
    const double feq_opp =
        vs_wi<VelocitySet>[iopp] * rho *
        (1.0 + 3.0 * cidotu_opp + 4.5 * cidotu_opp * cidotu_opp + omusq);

    // CALCULATE SYMMETRIC AND ANTISYMMETRIC PARTS OF THE DISTRIBUTION
    // FUNCTION
    const double fplus = 0.5 * (f[i] + f[iopp]);
    const double fminus = 0.5 * (f[i] - f[iopp]);
    const double fplus_eq = 0.5 * (feq_i + feq_opp);
    const double fminus_eq = 0.5 * (feq_i - feq_opp);

    // RELAX TO EQUILIBRIUM
    f[i] = f[i] - params.s_plus * (fplus - fplus_eq) -
           params.s_minus * (fminus - fminus_eq);
    f[iopp] = f[iopp] - params.s_plus * (fplus - fplus_eq) +
              params.s_minus * (fminus - fminus_eq);
  }
}

template <types::dim_t dim, typename VelocitySet, CollisionModel cm_t>
static __device__ __forceinline__ void
collide_node(double *__restrict__ f, const types::Coordinate<dim> p,
             const double rho, const utils::Vector<double, dim> u,
             const Params<dim, cm_t> params) {
  if constexpr (cm_t == lbm::CollisionModel::BGK) {
    collide_bgk_node(f, p, rho, u, params);
  } else if constexpr (cm_t == lbm::CollisionModel::TRT) {
    collide_trt_node(f, p, rho, u, params);
  } else {
    static_assert(assertion::always_false<dim>,
                  "collide_node: modello non implementato per CUDA.");
  }
}

} // namespace cuda
} // namespace lbm

#endif // __LBM_SIM_COLLISION_OPERATORS_COLLISION_STRATEGY_CUDA_CUH
