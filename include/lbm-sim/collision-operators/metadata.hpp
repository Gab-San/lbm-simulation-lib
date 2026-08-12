#ifndef __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP
#define __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

namespace lbm {
enum CollisionModel { BGK, TRT, MRT };

template <unsigned short int dim, enum CollisionModel cm_t> struct Params;
// constexpr bool always_false = false;
template <unsigned short int dim> struct Params<dim, CollisionModel::BGK> {
  const utils::Vector<double, dim> init_vel;
  const types::DimPoint<dim> num_cells;

  const double reyn_num;
  const double nu;
  utils::Vector<double, dim> F{6.66e-7, 0.0}; // default: nessuna forzante
  const double tauinv, omtauinv;

  Params(const double reyn_num_, const types::DimPoint<dim> num_cells_,
         const utils::Vector<double, dim> init_vel_)
      : init_vel(init_vel_), num_cells(num_cells_), reyn_num(reyn_num_),
        nu([&]() -> double {
          // WARN: this might need correction for different possible
          // configurations of velocity and cells
          if constexpr (dim == 2) {
            return init_vel_.dx * num_cells_.y / reyn_num_;
          } else {
            static_assert(dim != dim, "BGK : 3D not implemented yet!");
            // throw std::runtime_error("BGK: 3D not implemented yet!");
          }
        }()),
        tauinv(2.0 / (6.0 * nu + 1.0)), omtauinv(1.0 - tauinv) {
    double tau = 0.5 + 3.0 * nu;
    if (tau <= 0.5)
      throw std::runtime_error("LBM error: tau must be > 0.5");

    if (tau < 0.55 || tau > 1.2)
      std::cerr << "LBM warning: tau out of stability range, simulation may be "
                   "unstable."
                << std::endl;
  }
};

template <unsigned short int dim> struct Params<dim, CollisionModel::TRT> {
  const utils::Vector<double, dim> init_vel;
  const types::DimPoint<dim> num_cells;

  const double reyn_num;
  const double nu;

  const double tauPlus, tauMinus;
  const double s_plus, s_minus;

  Params(const double reyn_num_, const types::DimPoint<dim> num_cells_,
         const utils::Vector<double, dim> init_vel_)
      : init_vel(init_vel_), num_cells(num_cells_), reyn_num(reyn_num_),
        nu([&]() -> double {
          // WARN: this might need correction for different possible
          // configurations of velocity and cells
          if constexpr (dim == 2) {
            return init_vel_.dx * num_cells_.y / reyn_num_;
          } else {
            static_assert(dim != dim, "TRT : 3D not implemented yet!");
            // throw std::runtime_error("TRT: 3D not implemented yet!");
          }
        }()),
        tauPlus(3.0 * nu + 0.5), tauMinus(0.5 + (0.25) / (tauPlus - 0.5)),
        s_plus(1.0 / tauPlus), s_minus(1.0 / tauMinus) {
#if DEBUG
    // --- AGGIUNGI QUESTA RIGA ---
    std::cerr << "DEBUG: tauPlus=" << tauPlus << " tauMinus=" << tauMinus
              << " s_plus=" << s_plus << " s_minus=" << s_minus << std::endl;
    // -----------------------------
#endif

    if (tauPlus <= 0.5) {
      throw std::runtime_error("LBM error: tau must be > 0.5");
    }
  }
};

} // namespace lbm

#endif // __LBM_SIM_COLLISION_OPERATORS_METADATA_HPP
