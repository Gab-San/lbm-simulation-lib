/**
 * @file constants.hpp
 * @brief Numerical constants fixed by the LBM discretisation.
 *
 * Header-only, dependency-free. The values here are properties of the
 * lattice itself, not of a particular simulation, and are therefore
 * shared by every velocity set, collision model and backend.
 */

#ifndef __LBM_SIM_CONSTANTS_HPP
#define __LBM_SIM_CONSTANTS_HPP

namespace lbm {

/**
 * @brief Compile-time constants of the lattice discretisation.
 */
namespace numbers {

/**
 * @brief Inverse squared lattice speed of sound: @f$1 / c_s^2@f$.
 *
 * For the standard DdQq velocity sets with unit lattice spacing and unit
 * time step, @f$c_s^2 = 1/3@f$, hence the value 3.0.
 *
 * The constant appears in the equilibrium distribution and in the moment
 * computations:
 * @f[
 *   f_i^{eq} = w_i \rho \left( 1
 *     + \frac{c_i \cdot u}{c_s^2}
 *     + \frac{(c_i \cdot u)^2}{2 c_s^4}
 *     - \frac{u \cdot u}{2 c_s^2} \right)
 * @f]
 *
 * @note Declared @c static at namespace scope, so each translation unit
 *       gets its own internal-linkage copy; being @c constexpr, the value
 *       is normally folded away entirely.
 */
static constexpr double invcs_2 = 3.0;
} // namespace numbers
} // namespace lbm

#endif // __LBM_SIM_CONSTANTS_HPP
