/**
 * @file collision-params.hpp
 * @brief CollisionParams<dim, cm_t>: the relaxation constants of one
 *        collision operator, derived once from the physical inputs.
 *
 * The primary template is only declared. Each collision model provides its
 * own specialisation with the rates *that model* needs, so a parameter set
 * can never be paired with the wrong operator: the mismatch is a
 * substitution failure, not a run-time surprise.
 *
 * Every specialisation takes the same three physical inputs -- Reynolds
 * number, grid size and reference velocity -- and derives the kinematic
 * viscosity the same way:
 *
 * @f[
 *   \nu = \frac{u_0 \, N_y}{\mathrm{Re}}
 * @f]
 *
 * The characteristic length is @c num_cells.y in both 2D and 3D (the cavity
 * height, the channel height), and the characteristic velocity is
 * @c init_vel.dx. That is why @c init_vel.dx must be strictly positive even
 * in problems where no wall moves -- see the "Configuration files" page.
 *
 * @note The constructors are @c LBM_HD_FUNC so a parameter object can be
 *       built on the device, but the validation is fenced behind
 *       @c __CUDA_ARCH__: there is no throwing and no @c iostream in device
 *       code. Build the object on the host, where the checks run, and copy
 *       it down.
 */

#pragma once

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/constants.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/common.hpp"

#include <cstddef>

namespace lbm {

/**
 * @brief Relaxation parameters of a collision operator.
 *
 * Declared, never defined: only the specialisations below exist. Naming a
 * model without one -- @c CollisionModel::MRT today -- is an incomplete type
 * and fails to compile at the point of use.
 *
 * @tparam dim  Spatial dimension (2 or 3).
 * @tparam cm_t Collision model selecting the specialisation.
 */
template <unsigned short int dim, enum CollisionModel cm_t>
struct CollisionParams;

/**
 * @brief BGK parameters: a single relaxation time.
 *
 * @f[
 *   \tau = \tfrac{1}{2} + 3\nu, \qquad
 *   \omega = \frac{1}{\tau} = \frac{2}{6\nu + 1}
 * @f]
 *
 * @c tauinv and @c omtauinv are stored rather than @c tau because the
 * collision step wants @f$ f \leftarrow (1-\omega) f + \omega f^{eq} @f$:
 * keeping both forms turns the inner loop into two multiplies and an add,
 * with no division.
 */
template <unsigned short int dim>
struct CollisionParams<dim, CollisionModel::BGK> {
  /// Reference velocity. @c dx is the characteristic velocity of the problem.
  const utils::Vector<double, dim> init_vel;

  /// Reynolds number the viscosity was derived from.
  const double reyn_num;

  /// Kinematic viscosity in lattice units, @f$ \nu = u_0 N_y / Re @f$.
  const double nu;

  /// @f$ \omega = 1/\tau @f$ and @f$ 1 - \omega @f$, the two factors of the
  /// BGK relaxation.
  const double tauinv, omtauinv;

  /**
   * @brief Derives @c nu and the relaxation rates from the physical inputs.
   *
   * @param reyn_num_  Reynolds number, must be > 0.
   * @param num_cells_ Grid extents.
   * @param init_vel_  Reference velocity; @c dx must be > 0.
   *
   * @throws std::runtime_error if the resulting @f$ \tau \le 0.5 @f$, where
   *         the scheme has no positive viscosity and cannot be stable.
   *
   * @warning Prints a warning on @c stderr when @f$ \tau @f$ falls outside
   *          the roughly @c 0.55 -- @c 1.2 band. Outside it the simulation
   *          often still runs, and often still diverges.
   *
   * @note Both checks are skipped in device code (@c __CUDA_ARCH__).
   */
  LBM_HD_FUNC
  CollisionParams(const double reyn_num_,
                  const utils::Vector<double, dim> init_vel_,
                  const double u_ref, const std::size_t L)
      : init_vel(init_vel_), reyn_num(reyn_num_), nu(u_ref * L / reyn_num_),
        tauinv(2.0 / (2 * numbers::invcs_2 * nu + 1.0)),
        omtauinv(1.0 - tauinv) {
#ifndef __CUDA_ARCH__
    double tau = 0.5 + numbers::invcs_2 * nu;
    if (tau <= 0.5) {
      throw std::runtime_error("LBM error: tau must be > 0.5");
    }

    if (tau < 0.55 || tau > 1.2) {
      std::cerr << "LBM warning: tau out of stability range, simulation may be "
                   "unstable."
                << std::endl;
    }
#endif
  }
};

/**
 * @brief TRT parameters: separate rates for the symmetric and antisymmetric
 *        parts of the populations.
 *
 * @f$ \tau^+ @f$ carries the viscosity exactly as in BGK, while
 * @f$ \tau^- @f$ is fixed by the "magic number"
 * @f$ \Lambda = (\tau^+ - \tfrac12)(\tau^- - \tfrac12) @f$:
 *
 * @f[
 *   \tau^+ = \tfrac{1}{2} + 3\nu, \qquad
 *   \tau^- = \tfrac{1}{2} + \frac{\Lambda}{\tau^+ - \tfrac{1}{2}},
 *   \qquad \Lambda = \tfrac{1}{4}
 * @f]
 *
 * @f$ \Lambda = 1/4 @f$ is the classic choice: it places a halfway
 * bounce-back wall exactly midway between the last fluid node and the first
 * solid node independently of @f$ \nu @f$, which is what removes the
 * viscosity-dependent wall offset BGK suffers from. It is hard-coded here;
 * making it configurable would mean adding it to the constructor.
 */
template <unsigned short int dim>
struct CollisionParams<dim, CollisionModel::TRT> {
  /// Reference velocity. @c dx is the characteristic velocity of the problem.
  const utils::Vector<double, dim> init_vel;

  /// Reynolds number the viscosity was derived from.
  const double reyn_num;

  /// Kinematic viscosity in lattice units, @f$ \nu = u_0 N_y / Re @f$.
  const double nu;

  /// Relaxation times of the symmetric and antisymmetric parts.
  const double tauPlus, tauMinus;

  /// Their reciprocals, the rates the collision kernel actually applies.
  const double inv_plus, inv_minus;

  /**
   * @brief Derives @c nu, the two relaxation times and their rates.
   *
   * @param reyn_num_  Reynolds number, must be > 0.
   * @param num_cells_ Grid extents.
   * @param init_vel_  Reference velocity; @c dx must be > 0.
   *
   * @throws std::runtime_error if @f$ \tau^+ \le 0.5 @f$.
   *
   * @note Unlike the BGK specialisation, no warning is emitted for a
   *       @f$ \tau^+ @f$ outside the comfortable band: TRT is stable over a
   *       distinctly wider range, which is one of the reasons to pick it.
   */
  LBM_HD_FUNC CollisionParams(const double reyn_num_,
                              const utils::Vector<double, dim> init_vel_,
                              const double u_ref, const std::size_t L,
                              const double lambda = 0.25)
      : init_vel(init_vel_), reyn_num(reyn_num_), nu(u_ref * L / reyn_num_),
        tauPlus(numbers::invcs_2 * nu + 0.5),
        tauMinus(0.5 + lambda / (tauPlus - 0.5)), inv_plus(1.0 / tauPlus),
        inv_minus(1.0 / tauMinus) {
#ifndef __CUDA_ARCH__
    if (tauPlus <= 0.5) {
      throw std::runtime_error("LBM error: tau must be > 0.5");
    }
#endif
  }
};

} // namespace lbm
