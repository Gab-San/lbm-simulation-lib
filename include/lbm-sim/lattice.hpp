/**
 * @file lattice.hpp
 * @brief Lattice<dim>: geometry, macroscopic fields and boundary
 *        description of a simulation.
 */

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

/**
 * @brief Full state of a simulation: geometry, macroscopic fields,
 *        obstacles and pressure boundary values.
 *
 * The distribution functions are deliberately *not* stored here: they
 * belong to the solver, which owns the populations and chooses their
 * layout (AoS/SoA, host/device) according to its backend. The lattice
 * holds what every backend and every post-processing step needs to agree
 * on.
 *
 * @tparam dim Spatial dimension (2 or 3).
 *
 * @note All members except @ref u and @ref rho are @c const, so the
 *       implicit copy- and move-assignment operators are deleted: a
 *       Lattice can be copied or moved into a new object but never
 *       reassigned in place, and cannot be stored in containers that
 *       require assignability.
 */
template <unsigned short int dim> struct Lattice {
  /// Domain geometry; @c getArea() gives the cell count, @c size the extents.
  const Grid<dim> grid;

  /// Macroscopic velocity per cell (@c grid.getArea() entries), written
  /// back by the solver at the end of a run.
  std::vector<utils::Vector<double, dim>> u;

  /// Reference density, fixed at 1.0 by the constructor.
  const double rho0;

  /// Macroscopic density per cell, initialised to @ref rho0 everywhere and
  /// written back by the solver.
  std::vector<double> rho;

  /// Per-cell obstacle id: @c types::FLUID marks fluid cells, any other
  /// value indexes @ref obstacles.
  const types::solid_mask_t solid_mask; // obstacle ids

  /// Obstacle table: id -> {boundary condition, wall velocity}.
  const std::vector<Solid::ObstacleData<dim>> obstacles; // id -> {bc, u_wall}

  /// Boundary condition applied on each side of the domain.
  const Solid::DomainBC<dim> domain_bc;

  /// Inlet pressure, for pressure-driven setups.
  const double pin;

  /// Outlet pressure, for pressure-driven setups.
  const double pout;

  /**
   * @brief Builds the grid and allocates the macroscopic fields.
   *
   * @ref u and @ref rho are sized to @c grid.getArea(); @ref rho0 is set to
   * 1.0 and @ref rho filled with it.
   *
   * @param grid_dim_    Domain extents.
   * @param solid_mask_  Per-cell obstacle ids (moved in).
   * @param obstacles_   Obstacle table indexed by id (moved in).
   * @param domain_bc_   Boundary condition of each domain side.
   * @param pin_         Inlet pressure. Defaults to 0.
   * @param pout_        Outlet pressure. Defaults to 0.
   *
   * @note `obstacles` may be empty: the pointer is only ever dereferenced
   *       when oid != types::FLUID, which cannot happen with an all-FLUID
   *       mask. An obstacle-free simulation therefore passes an empty
   *       vector, not a dummy entry.
   *
   * @warning The mask and the obstacle table are not cross-checked here.
   *          An id in @p solid_mask_ that is out of range for
   *          @p obstacles_ is undefined behaviour at solve time, not a
   *          construction-time error.
   */
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
