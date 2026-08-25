#ifndef __LBM_SIM_FUNCTIONS_HPP
#define __LBM_SIM_FUNCTIONS_HPP

#include "lbm-sim/lattice.hpp"

#include <vector>

namespace lbm {
namespace functional {

inline std::vector<double>
extract_dx_profile_along_y_center(const Lattice<2> &lattice) {
  std::vector<double> profile(lattice.grid.size.y);
  int x = lattice.grid.size.x / 2;
  for (auto y = 0; y < lattice.grid.size.y; ++y) {
    profile[y] = lattice.u[lattice.grid.scalar_index({x, y})].dx;
  }
  return profile;
}

inline std::vector<double>
extract_dy_profile_along_x_center(const Lattice<2> &lattice) {
  std::vector<double> profile(lattice.grid.size.x);
  int y = lattice.grid.size.y / 2;
  for (auto x = 0; x < lattice.grid.size.x; ++x) {
    profile[x] = lattice.u[lattice.grid.scalar_index({x, y})].dy;
  }
  return profile;
}

} // namespace functional
} // namespace lbm

#endif // __LBM_SIM_FUNCTIONS_HPP
