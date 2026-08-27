#ifndef _LBM_SIM_CORE_COLLISION_OPERATORS_HPP
#define _LBM_SIM_CORE_COLLISION_OPERATORS_HPP

#include "lbm-sim/annotations.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/operators.hpp"

#include "lbm-sim/backend.hpp"

namespace lbm {

template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
class CollisionStrategy {
public:
  const CollisionParams<dim, cm_t> params;
  static constexpr CollisionModel type = cm_t;

  LBM_HD_FUNC CollisionStrategy(const CollisionParams<dim, cm_t> &params)
      : params(params) {}

  LBM_HD_FUNC inline void apply(double *RESTRICT fp,
                                const types::Coordinate<dim> p,
                                const double localrho,
                                const utils::Vector<double, dim> u) const {
    if constexpr (cm_t == CollisionModel::BGK) {
      apply_bgk(fp, p, u, localrho);
    } else if constexpr (cm_t == CollisionModel::TRT) {
      apply_trt(fp, p, u, localrho);
    } else {
      static_assert(
          cm_t == CollisionModel::MRT,
          "CollisionStrategy::apply(), with cm_t = MRT not yet implemented!");
    }
  }

private:
  LBM_HD_FUNC inline void apply_bgk(double *RESTRICT fp,
                                    const types::Coordinate<dim> p,
                                    const utils::Vector<double, dim> u,
                                    const double localrho) const {
    using utils::ops::dot;
    const double omusq = -1.5 * dot(u, u);

#ifndef __CUDA_ARCH__
#pragma omp simd
#endif
    for (auto diridx = 0; diridx < VelocitySet::ndir; ++diridx) {
      const double cidotu = dot(detail::direction<dim, VelocitySet>(diridx), u);
      const double feq = detail::weight<VelocitySet>(diridx) * localrho *
                         (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

      // RELAX TO EQUILIBRIUM
      fp[diridx] = params.omtauinv * fp[diridx] + params.tauinv * feq;
    }
  }

  LBM_HD_FUNC inline void apply_trt(double *RESTRICT fp,
                                    const types::Coordinate<dim> p,
                                    const utils::Vector<double, dim> u,
                                    const double localrho) const {
    using utils::ops::dot;
    const double omusq = -1.5 * dot(u, u);

#ifndef __CUDA_ARCH__
#pragma omp simd
#endif
    for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
      const auto iopp = detail::opposite<VelocitySet>(i);

      // NOTE: WHY THIS CHECK?
      if (i > iopp) {
        continue;
      }

      const double cidotu_i = dot(detail::direction<dim, VelocitySet>(i), u);
      const double feq_i =
          detail::weight<VelocitySet>(i) * localrho *
          (1.0 + 3.0 * cidotu_i + 4.5 * cidotu_i * cidotu_i + omusq);

      if (i == iopp) {
        // TODO: THIS COMMENT IS SHIT ENGLISH
        //
        // Center Direction: no antisymmetric component.
        fp[i] = fp[i] - params.s_plus * (fp[i] - feq_i);
        continue;
      }

      const double cidotu_opp =
          dot(detail::direction<dim, VelocitySet>(iopp), u);

      // calculate equilibrium
      const double feq_opp =
          detail::weight<VelocitySet>(iopp) * localrho *
          (1.0 + 3.0 * cidotu_opp + 4.5 * cidotu_opp * cidotu_opp + omusq);

      // CALCULATE SYMMETRIC AND ANTISYMMETRIC PARTS OF THE DISTRIBUTION
      // FUNCTION
      const double fplus = 0.5 * (fp[i] + fp[iopp]);
      const double fminus = 0.5 * (fp[i] - fp[iopp]);
      const double fplus_eq = 0.5 * (feq_i + feq_opp);
      const double fminus_eq = 0.5 * (feq_i - feq_opp);

      // RELAX TO EQUILIBRIUM
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
