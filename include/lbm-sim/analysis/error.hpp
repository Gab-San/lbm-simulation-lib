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

/// CSV schema for an error table, in the shape CsvWriter expects.
/// @see the "Output formats" page.
struct ErrorAnalysisSchema {
  static constexpr char const *header =
      "profile,size,collision_model,niters,type,rel,abs";
  static constexpr char const *format = "{},{},{},{},{},{:3f},{:3f}";
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
template <unsigned short int dim> class ErrorEvaluator {
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

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_ERROR_HPP
