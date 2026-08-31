/**
 * @file collision-strategy.hpp
 * @brief CollisionStrategy: the collision kernel itself, selected at compile
 *        time by the collision model.
 *
 * One object, one @c apply() entry point, and an @c if @c constexpr chain
 * that resolves to a single operator at compile time. There is no virtual
 * dispatch and no branch on the model inside the loop over the directions:
 * the strategy is constructed once per solve and inlined into the fused
 * stream-collide pass.
 *
 * @see collision-params.hpp for the relaxation constants each operator uses.
 */

#ifndef _LBM_SIM_CORE_COLLISION_OPERATORS_HPP
#define _LBM_SIM_CORE_COLLISION_OPERATORS_HPP

#include "lbm-sim/annotations.hpp"
#include "lbm-sim/backend/utils.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"

#include "lbm-sim/core/operators.hpp"

namespace lbm {

/**
 * @brief Relaxes the populations of one node towards equilibrium.
 *
 * @tparam dim         Spatial dimension (2 or 3).
 * @tparam VelocitySet Discrete velocity set supplying @c ndir, @c dir,
 *                     @c wi and @c opp.
 * @tparam cm_t        Collision model; also selects the CollisionParams
 *                     specialisation held by @ref params.
 *
 * All members are @c LBM_HD_FUNC, so the same code compiles for the OpenMP
 * and the CUDA backends. Adding an operator therefore means writing device-
 * compatible code: no allocation, no exceptions, no standard containers.
 */
template <types::dim_t dim, typename VelocitySet, enum CollisionModel cm_t>
class CollisionStrategy {
public:
  /// Relaxation constants, copied in at construction.
  const CollisionParams<dim, cm_t> params;

  /// The model this strategy implements, readable from generic code.
  static constexpr CollisionModel type = cm_t;

  /**
   * @brief Stores a copy of the relaxation parameters.
   * @param params Parameter set matching @p cm_t.
   */
  LBM_HD_FUNC CollisionStrategy(const CollisionParams<dim, cm_t> &params)
      : params(params) {}

  /**
   * @brief Applies the collision operator in place.
   *
   * @param[in,out] fp       The @c VelocitySet::ndir populations of the node,
   *                         overwritten with the post-collision values.
   * @param[in]     p        Node coordinate. Unused by the current operators,
   *                         kept in the signature for space-dependent ones
   *                         (forcing terms, local viscosity).
   * @param[in]     localrho Macroscopic density at the node.
   * @param[in]     u        Macroscopic velocity at the node.
   *
   * Dispatch is an @c if @c constexpr chain closed by a @c static_assert, so
   * a model that has parameters but no kernel -- @c CollisionModel::MRT --
   * fails to compile where it is used, instead of silently doing nothing.
   */
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
  /**
   * @brief Single-relaxation-time (BGK) collision.
   *
   * @f[
   *   f_i \leftarrow (1-\omega)\, f_i + \omega\, f_i^{eq}(\rho, \mathbf{u})
   * @f]
   *
   * with the second-order equilibrium
   *
   * @f[
   *   f_i^{eq} = w_i \rho \left( 1 + 3\,\mathbf{c}_i\!\cdot\!\mathbf{u}
   *              + \tfrac{9}{2} (\mathbf{c}_i\!\cdot\!\mathbf{u})^2
   *              - \tfrac{3}{2} \mathbf{u}\!\cdot\!\mathbf{u} \right)
   * @f]
   *
   * One pass over the directions, no cross-direction dependency, so the loop
   * carries an @c omp @c simd on the host path.
   */
  LBM_HD_FUNC inline void apply_bgk(double *RESTRICT fp,
                                    const types::Coordinate<dim> p,
                                    const utils::Vector<double, dim> u,
                                    const double localrho) const {
    using utils::ops::dot;
    const double omusq = -1.5 * dot(u, u);

#ifndef __CUDA_ARCH__
#pragma omp simd
#endif
    for (std::ptrdiff_t diridx = 0;
         diridx < static_cast<std::ptrdiff_t>(VelocitySet::ndir); ++diridx) {
      const double cidotu = dot(detail::direction<dim, VelocitySet>(diridx), u);
      const double feq = detail::weight<VelocitySet>(diridx) * localrho *
                         (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu + omusq);

      // RELAX TO EQUILIBRIUM
      fp[diridx] = params.omtauinv * fp[diridx] + params.tauinv * feq;
    }
  }

  /**
   * @brief Two-relaxation-time (TRT) collision.
   *
   * Each population is split into a symmetric and an antisymmetric part with
   * respect to its opposite direction, and the two are relaxed at different
   * rates:
   *
   * @f[
   *   f_i^{\pm} = \tfrac{1}{2}\left( f_i \pm f_{\bar\imath} \right), \qquad
   *   f_i \leftarrow f_i - s^{+}\!\left(f_i^{+} - f_i^{eq,+}\right)
   *                      - s^{-}\!\left(f_i^{-} - f_i^{eq,-}\right)
   * @f]
   *
   * and symmetrically for @f$ f_{\bar\imath} @f$, which is why the loop
   * updates both members of a pair at once and skips the iteration when
   * @c i > @c iopp: each pair is visited exactly once, from its lower index.
   * The rest direction, where @c i == @c iopp, has no antisymmetric part and
   * reduces to a plain BGK step at rate @f$ s^{+} @f$.
   *
   * @warning Not yet validated against a reference to the same extent as the
   *          BGK path, and not optimised: the equilibrium of both directions
   *          of a pair is recomputed rather than shared, and the early
   *          @c continue keeps the loop from vectorising as cleanly as
   *          apply_bgk() does.
   *
   * @see CollisionParams<dim, CollisionModel::TRT> for how @f$ \tau^{\pm} @f$
   *      are derived and why @f$ \Lambda = 1/4 @f$.
   */
  LBM_HD_FUNC inline void apply_trt(double *RESTRICT fp,
                                    const types::Coordinate<dim> p,
                                    const utils::Vector<double, dim> u,
                                    const double localrho) const {
    using utils::ops::dot;
    const double omusq = -1.5 * dot(u, u);

    // The rest direction is its own opposite: no antisymmetric part, so
    // this degenerates to a BGK relaxation at rate inv_plus.
    const double cidotu_0 = dot(detail::direction<dim, VelocitySet>(0), u);
    const double feq_0 =
        detail::weight<VelocitySet>(0) * localrho *
        (1.0 + numbers::invcs_2 * cidotu_0 + 4.5 * cidotu_0 * cidotu_0 + omusq);
    fp[0] = fp[0] - params.inv_plus * (fp[0] - feq_0);

#ifndef __CUDA_ARCH__
#pragma omp simd
#endif
    for (std::ptrdiff_t i = 1;
         i < static_cast<std::ptrdiff_t>(VelocitySet::ndir); i += 2) {
      const auto iopp =
          static_cast<std::ptrdiff_t>(detail::opposite<VelocitySet>(i));

      const double cidotu_i = dot(detail::direction<dim, VelocitySet>(i), u);
      const double feq_i = detail::weight<VelocitySet>(i) * localrho *
                           (1.0 + numbers::invcs_2 * cidotu_i +
                            4.5 * cidotu_i * cidotu_i + omusq);

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
      fp[i] = fp[i] - params.inv_plus * (fplus - fplus_eq) -
              params.inv_minus * (fminus - fminus_eq);
      fp[iopp] = fp[iopp] - params.inv_plus * (fplus - fplus_eq) +
                 params.inv_minus * (fminus - fminus_eq);
    }
  }
};

} // namespace lbm

#endif
