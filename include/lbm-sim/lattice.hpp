#ifndef __LBM_SIM_LATTICE_HPP
#define __LBM_SIM_LATTICE_HPP

#include "lbm-sim/core/grid.hpp"

#include "collision-detection/core/vector.hpp"

// C++ STD LIB
#include <vector>

namespace lbm {
template <int dim> struct Lattice {

  const Grid<dim> grid;
  std::vector<CollisionDetection::utils::Vector<double, dim>> u;
  const double rho0;
  std::vector<double> rho;

  Lattice(CollisionDetection::types::DimPoint<dim> grid_dim_)
      : grid(grid_dim_), u(grid.getArea()), rho0(1.0),
        rho(grid.getArea(), rho0) {}
};
} // namespace lbm
#endif // __LBM_SIM_LATTICE_HPP
