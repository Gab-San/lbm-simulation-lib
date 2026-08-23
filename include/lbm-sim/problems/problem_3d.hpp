#pragma once

#include "lbm-sim/lattice.hpp"

#include "lbm-sim/boundaries.hpp"

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/operators.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <cstddef>
#include <vector>

#include <omp.h>

namespace lbm {

/**
 * Cavita' cubica 3D (analogo 3D di LidCavity2D, problems/problem_2d.hpp):
 * 5 pareti rigide (bounce-back semplice) + 1 parete mobile (il "lid",
 * bounce-back con velocita' imposta) sulla faccia superiore z = Nz-1.
 * La parete si muove lungo x, come nel caso 2D.
 *
 * NOTA: a differenza del caso 2D, qui non passiamo dal modulo
 * collision-detection per generare la boundary_mask: shape.hpp e'
 * ancora segnato "FIXME: Should be extendible to 3D" e Segment/Circle/
 * Parallelogram lavorano solo su (x,y). Per una cavita' cubica non
 * serve comunque quella generalita': la maschera si costruisce in modo
 * diretto scorrendo le facce del cubo (vedi build_boundary_mask sotto).
 */
class LidCavity3D {
public:
  virtual ~LidCavity3D() = default;

  /// Analogo a LidCavity2D::init(): imposta la velocita' iniziale sui
  /// nodi indicati (tipicamente i nodi del lid restituiti da
  /// lid_points()).
  void
  init(Lattice<3> &lattice, const utils::Vector<double, 3> &init_vel,
       const std::vector<types::Coordinate<3>> &initialization_points) const {
    for (const auto &p : initialization_points) {
      lattice.u[lattice.grid.scalar_index(p)] = init_vel;
    }
  }

  /**
   * Costruisce la boundary_mask per un dominio cubico Nx*Ny*Nz:
   *  - z == Nz-1                              -> BB_MOVING_WALL (il lid)
   *  - z == 0, x in {0,Nx-1}, y in {0,Ny-1}    -> BB_RIGID_WALL
   *  - nodi interni                           -> Solid::NONE
   */
  static types::boundary_mask_t
  build_boundary_mask(const types::DimPoint<3> &size) {
    types::boundary_mask_t mask(utils::ops::measure<3>(size), Solid::NONE);

#pragma omp parallel for collapse(3) schedule(static)
    for (std::size_t z = 0; z < size.z; ++z) {
      for (std::size_t y = 0; y < size.y; ++y) {
        for (std::size_t x = 0; x < size.x; ++x) {
          const bool on_boundary = (x == 0 || x == size.x - 1 || y == 0 ||
                                    y == size.y - 1 || z == 0 ||
                                    z == size.z - 1);
          if (!on_boundary)
            continue;

          const types::Coordinate<3> p(static_cast<int>(x),
                                       static_cast<int>(y),
                                       static_cast<int>(z));
          const std::size_t idx = size.x * (size.y * z + y) + x;

          mask[idx] = (z == size.z - 1) ? Solid::BB_MOVING_WALL
                                        : Solid::BB_RIGID_WALL;
        }
      }
    }

    return mask;
  }

  /// Nodi della faccia superiore (il lid) da passare a init(): sono
  /// esattamente i nodi marcati BB_MOVING_WALL da build_boundary_mask().
  static std::vector<types::Coordinate<3>>
  lid_points(const types::DimPoint<3> &size) {
    std::vector<types::Coordinate<3>> pts;
    pts.reserve(size.x * size.y);
    const std::size_t z = size.z - 1;
    for (std::size_t y = 0; y < size.y; ++y) {
      for (std::size_t x = 0; x < size.x; ++x) {
        pts.emplace_back(static_cast<int>(x), static_cast<int>(y),
                         static_cast<int>(z));
      }
    }
    return pts;
  }
};

} // namespace lbm
