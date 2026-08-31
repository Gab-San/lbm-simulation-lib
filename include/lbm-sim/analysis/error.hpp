/**
 * @file error.hpp
 * @brief ErrorEvaluator: the two-step error computation, in the style of
 *        dealii::VectorTools.
 *
 * The split is the same one deal.II makes, and for the same reason: a
 * per-cell error field first, then a single reduction over it. Keeping them
 * apart means the per-cell field can be inspected or exported (where is the
 * error, not just how big), and the norm can be changed without recomputing
 * the difference.
 *
 * @see LBMSimulation::compute_error(), which chains the two on the
 *      simulation's own lattice.
 */
#pragma once

#include "lbm-sim/analysis/types.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/types/base.hpp"

// C++ STANDARD LIB
#include <algorithm>
#include <cmath>
#include <vector>

namespace lbm {
namespace analysis {

/// CSV schema for an error table, in the shape CsvWriter expects.
/// @see the "Output formats" page.
struct ErrorAnalysisSchema {
  static constexpr char const *header =
      "profile,size,collision_model,niters,type,rel,abs";
  static constexpr char const *format = "{},{},{},{},{},{:3f},{:3f}";
};

/// CSV schema for an error table, in the shape CsvWriter expects.
/// @see the "Output formats" page.
struct DetailedErrorAnalysisSchema {
  static constexpr char const *header =
      "profile,grid_size,reynolds,collision_model,niters,type,rel,abs,"
      "rel_percent,rmse,rmse_ref_percent,linf,linf_ref_percent";
  static constexpr char const *format =
      "{},{},{},{},{},{},{:.8e},{:.8e},{:.6f},{:.8e},{:.6f},{:.8e},"
      "{:.6f}";
};

/**
 * @brief The two halves of an error computation: a per-cell field, then a
 *        reduction over it.
 *
 * @tparam dim Spatial dimension (2 or 3).
 *
 * Stateless; both members are static. LBMSimulation::compute_error() chains
 * them on its own lattice, which is the usual entry point.
 */
template <types::dim_t dim> class ErrorEvaluator {
public:
  /**
   * @brief Per-cell magnitude of the difference between the simulated and
   *        the exact velocity.
   *
   * @param grid           Domain, iterated in full.
   * @param sim_u          Simulated velocity field, @c grid.getArea()
   *                       entries.
   * @param exact_solution Analytical solution, sampled at each node's
   *                       integer coordinate.
   * @return One non-negative value per cell, unnormalised, indexed by
   *         Grid::scalar_index().
   *
   * @warning Every cell contributes, solid ones included. Where the solver
   *          leaves @c u = 0 inside a wall, the exact solution must do the
   *          same or the error is dominated by the interior of the obstacle
   *          -- which is why HagenPoiseuilleSolution3D returns zero outside
   *          the pipe instead of continuing the parabola.
   */
  static std::vector<double>
  integrate_difference(const Grid<dim> &grid,
                       const std::vector<utils::Vector<double, dim>> &sim_u,
                       const functional::Function<dim> &exact_solution) {

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
      // Same computation in 3D: the vector difference gains one component,
      // the rest (one unnormalised error per cell) is identical.
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

  /**
   * @brief Reduces a per-cell error field to a single number.
   *
   * @tparam T             Any range of @c double: the @c std::vector from
   *                       integrate_difference(), or the fixed-size array of
   *                       tabulated points the Ghia comparison builds.
   * @param error_per_cell Non-negative per-cell magnitudes.
   * @param norm_type      @c L1 (sum), @c L2 (root of the sum of squares),
   *                       @c L2_squared, or @c Linfty (maximum).
   * @return The norm, unnormalised.
   *
   * @note No cell-volume weighting: in lattice units the cell measure is 1,
   *       so an @c L1 norm is a plain sum. Comparing two different grid
   *       resolutions therefore requires normalising by the cell count.
   */
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
              const functional::Function<dim> &exact_solution,
              const NormType norm_type, const double reference_velocity = 0.0) {
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

  const double l2 =
      ErrorEvaluator<dim>::compute_global_error(error_per_cell, NormType::L2);
  const double linf = ErrorEvaluator<dim>::compute_global_error(
      error_per_cell, NormType::Linfty);
  const double rmse =
      error_per_cell.empty()
          ? 0.0
          : l2 / std::sqrt(static_cast<double>(error_per_cell.size()));

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
