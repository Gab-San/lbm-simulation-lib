/**
 * @file functions.hpp
 * @brief Ready-made profile extractors for LBMSimulation::output().
 *
 * Each reduces a lattice to the 1D velocity profile along a centreline --
 * the quantity the analytical solutions and the Ghia tables are compared
 * against. Any callable with the same signature works, so these are a
 * convenience, not an interface.
 *
 * @see the "Output formats" page for what output() does with the result.
 */

#ifndef __LBM_SIM_FUNCTIONS_HPP
#define __LBM_SIM_FUNCTIONS_HPP

#include "lbm-sim/lattice.hpp"

#include <cstddef>
#include <string_view>
#include <vector>

namespace lbm {
namespace functional {

template <types::dim_t dim> class Function {
public:
  virtual ~Function() = default;

  // Returns the exact velocity vector at the continuous point p.
  virtual utils::Vector<double, dim>
  value(const types::Coordinate<dim> &p) const = 0;

  virtual const std::string_view &getName() const = 0;
};

template <types::dim_t dim>
inline std::vector<double>
extract_dx_profile_along_y_center(const Lattice<dim> &lattice) {
  std::vector<double> profile(lattice.grid.size.y);
  const int x = lattice.grid.size.x / 2;
  if constexpr (dim == 2) {
    for (int y = 0; y < static_cast<int>(lattice.grid.size.y); ++y) {
      profile[y] = lattice.u[lattice.grid.scalar_index({x, y})].dx;
    }
  } else {
    const int z = lattice.grid.size.z / 2;
    for (int y = 0; y < static_cast<int>(lattice.grid.size.y); ++y) {
      profile[y] = lattice.u[lattice.grid.scalar_index({x, y, z})].dx;
    }
  }

  return profile;
}

template <types::dim_t dim>
inline std::vector<double>
extract_dy_profile_along_x_center(const Lattice<dim> &lattice) {
  std::vector<double> profile(lattice.grid.size.x);
  const int y = lattice.grid.size.y / 2;
  if constexpr (dim == 2) {
    for (int x = 0; x < static_cast<int>(lattice.grid.size.x); ++x) {
      profile[x] = lattice.u[lattice.grid.scalar_index({x, y})].dy;
    }
  } else {
    const int z = lattice.grid.size.z / 2;
    for (int x = 0; x < static_cast<int>(lattice.grid.size.x); ++x) {
      profile[x] = lattice.u[lattice.grid.scalar_index({x, y, z})].dy;
    }
  }
  return profile;
}

/**
 * \brief 3D counterpart of extract_dx_profile_along_y_center(): extracts ux
 * along the central vertical line of the domain (x = nx/2, y = ny/2), as z
 * varies.
 *
 * This is the profile the cubic lid cavity is usually compared against: the
 * lid moves along x on the z = nz-1 face, so ux(z) on the central column is
 * the direct counterpart of the 2D profile along y.
 */
inline std::vector<double>
extract_dx_profile_along_z_center(const Lattice<3> &lattice) {
  std::vector<double> profile(lattice.grid.size.z);
  const int x = lattice.grid.size.x / 2;
  const int y = lattice.grid.size.y / 2;
  for (int z = 0; z < static_cast<int>(lattice.grid.size.z); ++z) {
    profile[z] = lattice.u[lattice.grid.scalar_index({x, y, z})].dx;
  }
  return profile;
}

} // namespace functional
} // namespace lbm

#endif // __LBM_SIM_FUNCTIONS_HPP
