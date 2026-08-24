#ifndef __LBM_SIM_ANALYSIS_ERROR_HPP
#define __LBM_SIM_ANALYSIS_ERROR_HPP

// LBM SIM LIB
#include "lbm-sim/core/grid.hpp"

#include "lbm-sim/analysis/types.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <cmath>
#include <vector>

namespace lbm {
namespace analysis {
template <unsigned short int dim> class ErrorEvaluator {
public:
  static std::vector<double>
  integrate_difference(const Grid<dim> &grid,
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
                    "ErrorEvaluator::integrate_difference() non ancora "
                    "implementato per dim > 2!");
    }

    return error_per_cell;
  }

  template <typename T>
  static double compute_global_error(const T &error_per_cell,
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
