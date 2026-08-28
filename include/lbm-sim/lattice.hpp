#ifndef __LBM_SIM_LATTICE_HPP
#define __LBM_SIM_LATTICE_HPP

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/vector.hpp"

// boundaries/types.hpp, NOT boundary-conditions.hpp: DomainBC and ObstacleData
// without dragging collision-area.hpp and <omp.h> into every translation unit.
#include "lbm-sim/boundaries/types.hpp"
#include "lbm-sim/types/common.hpp"

// C++ STD LIB
#include <utility>
#include <vector>

namespace lbm {
template <unsigned short int dim> struct Lattice {
  const Grid<dim> grid;

  std::vector<utils::Vector<double, dim>> u;
  const double rho0;
  std::vector<double> rho;

  const types::solid_mask_t solid_mask;                  // obstacle ids
  const std::vector<Solid::ObstacleData<dim>> obstacles; // id -> {bc, u_wall}
  const Solid::DomainBC<dim> domain_bc;

  const double pin;
  const double pout;

  /// `obstacles` may be empty: the pointer is only ever dereferenced when
  /// oid != types::FLUID, which cannot happen with an all-FLUID mask.
  Lattice(types::DimPoint<dim> grid_dim_, types::solid_mask_t solid_mask_,
          std::vector<Solid::ObstacleData<dim>> obstacles_,
          Solid::DomainBC<dim> domain_bc_, const double pin_ = 0,
          const double pout_ = 0)
      : grid(grid_dim_), u(grid.getArea()), rho0(1.0),
        rho(grid.getArea(), rho0), solid_mask(std::move(solid_mask_)),
        obstacles(std::move(obstacles_)), domain_bc(domain_bc_), pin(pin_),
        pout(pout_) {}
};
} // namespace lbm

#endif // __LBM_SIM_LATTICE_HPP
