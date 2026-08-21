#ifndef __LBM_SIM_ANALYSIS_ERROR_HPP
#define __LBM_SIM_ANALYSIS_ERROR_HPP

// LBM SIM LIB
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/lattice.hpp"

// COLLISION DETECTION LIB
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace lbm {
namespace analysis {

enum class NormType { L1, L2, L2_squared, Linfty };

inline const char *to_string(NormType t) {
  switch (t) {
  case NormType::L1:
    return "L1";
  case NormType::L2:
    return "L2";
  case NormType::L2_squared:
    return "L2^2";
  case NormType::Linfty:
    return "Linf";
  }
  return "unknown";
}


/**
 * Risultato di un confronto con una soluzione analitica generica
 * (Couette, Poiseuille, ...) dove la norma è scelta dal chiamante.
 */
struct NormErrorResult {
  double relative;
  double absolute;
  NormType norm_type;
};

template <unsigned short int dim> class Function {
public:
  virtual ~Function() = default;

  // Restituisce la velocità vettoriale esatta nel punto continuo P
  virtual utils::Vector<double, dim>
  value(const types::Coordinate<dim> &p) const = 0;
};

template <unsigned short int dim> class ErrorEvaluator {
public:
  static std::vector<double> integrate_difference(
      const Grid<dim> &grid,
      const std::vector<utils::Vector<double, dim>> &sim_u,
      const Function<dim> &exact_solution) {

    std::vector<double> error_per_cell(grid.getArea(), 0.0);

    if constexpr (dim == 2) {
      for (int y = 0; y < static_cast<int>(grid.size.y); ++y) {
        for (int x = 0; x < static_cast<int>(grid.size.x); ++x) {
          types::Coordinate<dim> coord{x, y};
          std::size_t idx = grid.scalar_index(coord);

          utils::Vector<double, dim> u_exact = exact_solution.value(coord);
          utils::Vector<double, dim> u_sim = sim_u[idx];

          utils::Vector<double, dim> diff = u_sim - u_exact;
          error_per_cell[idx] =
              std::sqrt(diff.dx * diff.dx + diff.dy * diff.dy);
        }
      }
    } else {
      static_assert(assertion::always_false<dim>,
                    "ErrorEvaluator::integrate_difference() non ancora implementato per dim > 2!");
    }

    return error_per_cell;
  }

  static double compute_global_error(const std::vector<double> &error_per_cell,
                                     NormType norm_type) {
    double total_error = 0.0;
    double max_error = 0.0;

    for (double err : error_per_cell) {
      switch (norm_type) {
      case NormType::L1:
        total_error += err;
        break;
      case NormType::L2:
      case NormType::L2_squared:
        total_error += err * err;
        break;
      case NormType::Linfty:
        max_error = std::max(max_error, err);
        break;
      }
    }

    switch (norm_type) {
    case NormType::L1:
      return total_error;
    case NormType::L2_squared:
      return total_error;
    case NormType::L2:
      return std::sqrt(total_error);
    case NormType::Linfty:
      return max_error;
    default:
      return 0.0;
    }
  }
};


} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_ERROR_HPP