#ifndef __LBM_SIM_LATTICE_HPP
#define __LBM_SIM_LATTICE_HPP

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STD LIB
#include <vector>

namespace lbm {
template <unsigned short int dim> struct Lattice {
  const Grid<dim> grid;

  std::vector<utils::Vector<double, dim>> u;
  const double rho0;
  std::vector<double> rho;

  const types::boundary_mask_t boundary_mask;
  const double pin;
  const double pout;

  Lattice(types::DimPoint<dim> grid_dim_, types::boundary_mask_t boundary_mask_,
          const double pin_ = 0, const double pout_ = 0)
      : grid(grid_dim_), u(grid.getArea()), rho0(1.0),
        rho(grid.getArea(), rho0), boundary_mask(std::move(boundary_mask_)),
        pin(pin_), pout(pout_) {}
};
} // namespace lbm

#endif // __LBM_SIM_LATTICE_HPP
