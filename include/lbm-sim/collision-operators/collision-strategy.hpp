#ifndef _LBM_SIM_CORE_COLLISION_OPERATORS_HPP
#define _LBM_SIM_CORE_COLLISION_OPERATORS_HPP

#include "lbm-sim/backend/metadata.hpp"

#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/types.hpp"

namespace lbm {

template <unsigned short int dim, typename VelocitySet,
          enum CollisionModel cm_t, enum ExecutionBackend backend_t>
class CollisionStrategy {
public:
  const Params<dim, cm_t> params;
  static constexpr CollisionModel type = cm_t;

  CollisionStrategy(const Params<dim, cm_t> &params) : params(params) {}

  void apply(const CollisionDetection::types::Coordinate<dim> p,
             const CollisionDetection::utils::Vector<double, dim> u,
             const double localrho, std::array<double, VelocitySet::ndir> &fp,
             const Grid<dim> &grid) const {
    if constexpr (cm_t == CollisionModel::BGK) {
      apply_bgk(p, u, localrho, fp);
    } else if constexpr (cm_t == CollisionModel::TRT) {
      apply_trt(p, u, localrho, fp);
    } else {
      static_assert(
          cm_t == CollisionModel::MRT,
          "CollisionStrategy::apply(), with cm_t = MRT not yet implemented!");
    }
  }

private:
  void apply_bgk(const CollisionDetection::types::Coordinate<dim> p,
                 const CollisionDetection::utils::Vector<double, dim> u,
                 const double localrho,
                 std::array<double, VelocitySet::ndir> &fp) const {
    using CollisionDetection::utils::dot;

    const double omusq = -1.5 * dot(u, u);
    const double force_prefactor = 1.0 - 0.5 * params.tauinv;
    // Collisione con SIMD
#pragma omp simd
    for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
      // calculate dot product beetwen the velocity u(x,y)
      // and the direction vector to its neighbour
      const double cidotu = dot(VelocitySet::dir[i], u);

      // calculate equilibrium
      const double feq = VelocitySet::wi[i] * localrho *
                         (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);
      // relax to equilibrium
      int pois = 0;
      if (pois == 1) {
        const double cidotF = dot(VelocitySet::dir[i], params.F);
        const double udotF = dot(u, params.F);
        const double force_term =
            force_prefactor * VelocitySet::wi[i] *
            (3.0 * (cidotF - udotF) + 9.0 * cidotu * cidotF);

        // ---------------------------------------

        fp[i] = params.omtauinv * fp[i] + params.tauinv * feq + force_term;
      } else {
        fp[i] = params.omtauinv * fp[i] + params.tauinv * feq;
      }
    }
  }

  void apply_trt(const CollisionDetection::types::Coordinate<dim> p,
                 const CollisionDetection::utils::Vector<double, dim> u,
                 const double localrho,
                 std::array<double, VelocitySet::ndir> &fp) const {
    using CollisionDetection::utils::dot;

    const double omusq = -1.5 * dot(u, u);

    // Collisione con SIMD
#pragma omp simd
    for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
      const auto iopp = VelocitySet::opp[i];

      // NOTE: WHY THIS CHECK?
      if (i > iopp) {
        continue;
      }

      const double cidotu_i = dot(VelocitySet::dir[i], u);
      const double feq_i =
          VelocitySet::wi[i] * localrho *
          (1.0 + 3.0 * cidotu_i + 4.5 * cidotu_i * cidotu_i + omusq);

      if (i == iopp) {
        // direzione di riposo: nessuna componente antisimmetrica, un solo
        // update

        fp[i] = fp[i] - params.s_plus * (fp[i] - feq_i);

        continue;
      }

      const double cidotu_opp = dot(VelocitySet::dir[iopp], u);

      // calculate equilibrium
      const double feq_opp =
          VelocitySet::wi[iopp] * localrho *
          (1.0 + 3.0 * cidotu_opp + 4.5 * cidotu_opp * cidotu_opp + omusq);

      // calculate symmetric and antisymmetric parts of the distribution
      // function
      const double fplus = 0.5 * (fp[i] + fp[iopp]);
      const double fminus = 0.5 * (fp[i] - fp[iopp]);
      const double fplus_eq = 0.5 * (feq_i + feq_opp);
      const double fminus_eq = 0.5 * (feq_i - feq_opp);

      // relax to equilibrium
      fp[i] = fp[i] - params.s_plus * (fplus - fplus_eq) -
              params.s_minus * (fminus - fminus_eq);
      fp[iopp] = fp[iopp] - params.s_plus * (fplus - fplus_eq) +
                 params.s_minus * (fminus - fminus_eq);

      /*SOME IMPORTANT REMARKS*/
      /*
      the TRT collision operator
      -me must choose a value for tau_plus, which is the relaxation time for the
      symmetric part of the distribution function. ù This value is typically
      chosen to be close to 1.0, which corresponds to a low viscosity fluid.
      -Once tau_plus is chosen, we can calculate tau_minus using the relation
      tau_minus = 0.5 + 0.25/(tau_plus - 0.5). This ensures that the relaxation
      times are consistent with the desired viscosity of the fluid. -we must
      discuss on the stability and the choosing of tau_plus in base of what we
      need -the code isn't tested yed and optimized, so it may be not correct
      and/or not efficient. -the for can be optimized by calculating the
      equilibrium only once for each direction and reusing it for both the
      symmetric and antisymmetric parts of the distribution function.
      */
    }
  }
};

} // namespace lbm

#endif
