#ifndef __LBM_SIM_ANALYSIS_ERROR_HPP
#define __LBM_SIM_ANALYSIS_ERROR_HPP

// LBM SIM LIB
#include "lbm-sim/core/grid.hpp"

#include "lbm-sim/analysis/types.hpp"
#include "lbm-sim/core/vector.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <cmath>
#include <vector>

namespace lbm {
namespace analysis {

struct ErrorAnalysisSchema {
  static constexpr char const *header =
      "profile,size,collision_model,niters,type,rel,abs";
  static constexpr char const *format = "{},{},{},{},{},{:3f},{:3f}";
};

//rmse, Linf, normalized percentages
// Grid dimensions are stored in separate columns: writing DimPoint directly
// would inject the comma from "Point(Nx,Ny)" into the CSV row.
struct DetailedErrorAnalysisSchema {
  static constexpr char const *header =
      "profile,grid_x,grid_y,reynolds,collision_model,niters,type,rel,abs,"
      "rel_percent,rmse,rmse_ref_percent,linf,linf_ref_percent,runtime_s,"
      "threads,mlups";
  static constexpr char const *format =
      "{},{},{},{},{},{},{},{:.8e},{:.8e},{:.6f},{:.8e},{:.6f},{:.8e},"
      "{:.6f},{:.6f},{},{:.6f}";
};

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
      // Stesso conto in 3D: la differenza vettoriale ha una componente in
      // piu', il resto (un errore per cella, non normalizzato) e' identico.
      for (int z = 0; z < static_cast<int>(grid.size.z); ++z) {
        for (int y = 0; y < static_cast<int>(grid.size.y); ++y) {
          for (int x = 0; x < static_cast<int>(grid.size.x); ++x) {
            types::Coordinate<dim> coord{x, y, z};
            std::size_t idx = grid.scalar_index(coord);

            utils::Vector<double, dim> u_exact = exact_solution.value(coord);
            utils::Vector<double, dim> u_sim = sim_u[idx];

            utils::Vector<double, dim> diff = u_sim - u_exact;
            error_per_cell[idx] = std::sqrt(
                diff.dx * diff.dx + diff.dy * diff.dy + diff.dz * diff.dz);
          }
        }
      }
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

template <unsigned short int dim>
inline NormErrorResult
compute_error(const Grid<dim> &grid,
              const std::vector<utils::Vector<double, dim>> &sim_u,
              const Function<dim> &exact_solution, const NormType norm_type,
              const double reference_velocity = 0.0) {
  const auto error_per_cell =
      ErrorEvaluator<dim>::integrate_difference(grid, sim_u, exact_solution);

  std::vector<double> reference_per_cell(grid.getArea(), 0.0);

  for (std::size_t idx = 0; idx < grid.getArea(); ++idx) {
    types::Coordinate<dim> coord;

    if constexpr (dim == 2) {
      coord.x = static_cast<int>(idx % grid.size.x);
      coord.y = static_cast<int>(idx / grid.size.x);
    } else {
      const std::size_t plane = grid.size.x * grid.size.y;
      coord.z = static_cast<int>(idx / plane);
      const std::size_t rem = idx % plane;
      coord.y = static_cast<int>(rem / grid.size.x);
      coord.x = static_cast<int>(rem % grid.size.x);
    }

    const auto u_exact = exact_solution.value(coord);
    if constexpr (dim == 2) {
      reference_per_cell[idx] =
          std::sqrt(u_exact.dx * u_exact.dx + u_exact.dy * u_exact.dy);
    } else {
      reference_per_cell[idx] =
          std::sqrt(u_exact.dx * u_exact.dx + u_exact.dy * u_exact.dy +
                    u_exact.dz * u_exact.dz);
    }
  }

  const double absolute =
      ErrorEvaluator<dim>::compute_global_error(error_per_cell, norm_type);
  const double reference_norm =
      ErrorEvaluator<dim>::compute_global_error(reference_per_cell, norm_type);

  const double relative =
      reference_norm > 0.0 ? absolute / reference_norm : 0.0;

  const double l2 = ErrorEvaluator<dim>::compute_global_error(
      error_per_cell, NormType::L2);
  const double linf = ErrorEvaluator<dim>::compute_global_error(
      error_per_cell, NormType::Linfty);
  const double rmse = error_per_cell.empty()
                          ? 0.0
                          : l2 / std::sqrt(
                                     static_cast<double>(error_per_cell.size()));

  const double u_ref = std::abs(reference_velocity);

  NormErrorResult result{relative, absolute, norm_type};
  result.rmse = rmse;
  result.linf = linf;
  result.rmse_normalized = u_ref > 0.0 ? rmse / u_ref : 0.0;
  result.linf_normalized = u_ref > 0.0 ? linf / u_ref : 0.0;
  result.sample_count = error_per_cell.size();
  return result;
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_ERROR_HPP
