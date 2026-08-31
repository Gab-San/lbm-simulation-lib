#ifndef __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP
#define __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP

#include "lbm-sim/backend/cuda/annotations.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/common.hpp"

namespace lbm {

template <unsigned short int dim, enum CollisionModel cm_t>
struct CollisionParams;

template <unsigned short int dim>
struct CollisionParams<dim, CollisionModel::BGK> {
  const utils::Vector<double, dim> init_vel;
  const types::DimPoint<dim> num_cells;

  const double reyn_num;
  const double nu;
  const double tauinv, omtauinv;

  LBM_HD_FUNC CollisionParams(const double reyn_num_,
                              const types::DimPoint<dim> num_cells_,
                              const utils::Vector<double, dim> init_vel_)
      : init_vel(init_vel_), num_cells(num_cells_), reyn_num(reyn_num_),
        // Caratteristica: velocità del lid (init_vel.dx) e altezza cavità
        // (num_cells.y), valide indipendentemente da dim (2D o 3D).
        nu(init_vel_.dx * num_cells_.y / reyn_num_),
        tauinv(2.0 / (6.0 * nu + 1.0)), omtauinv(1.0 - tauinv) {
#ifndef __CUDA_ARCH__
    double tau = 0.5 + 3.0 * nu;
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

template <unsigned short int dim>
struct CollisionParams<dim, CollisionModel::TRT> {
  const utils::Vector<double, dim> init_vel;
  const types::DimPoint<dim> num_cells;

  const double reyn_num;
  const double nu;

  const double tauPlus, tauMinus;
  const double s_plus, s_minus;

  LBM_HD_FUNC CollisionParams(const double reyn_num_,
                              const types::DimPoint<dim> num_cells_,
                              const utils::Vector<double, dim> init_vel_)
      : init_vel(init_vel_), num_cells(num_cells_), reyn_num(reyn_num_),
        nu(init_vel_.dx * num_cells_.y / reyn_num_), tauPlus(3.0 * nu + 0.5),
        tauMinus(0.5 + (0.25) / (tauPlus - 0.5)), s_plus(1.0 / tauPlus),
        s_minus(1.0 / tauMinus) {
#if DEBUG
    std::cerr << "DEBUG: tauPlus=" << tauPlus << " tauMinus=" << tauMinus
              << " s_plus=" << s_plus << " s_minus=" << s_minus << std::endl;
#endif

#ifndef __CUDA_ARCH__
    if (tauPlus <= 0.5) {
      throw std::runtime_error("LBM error: tau must be > 0.5");
    }
#endif
  }
};

} // namespace lbm

#endif // __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP
