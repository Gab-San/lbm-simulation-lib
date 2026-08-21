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

struct ErrorResult {
  double l2_relative;
  double linf_absolute;
};

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
  value(const types::ValuePoint<dim> &p) const = 0;
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

          types::ValuePoint<dim> phys_point{static_cast<double>(x),
                                            static_cast<double>(y)};

          utils::Vector<double, dim> u_exact = exact_solution.value(phys_point);
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

  // Come sopra, ma per la magnitudine della sola soluzione esatta:
  // serve come termine di normalizzazione in compute_error(), con
  // la stessa norma usata per l'errore (altrimenti il rapporto
  // relativo non ha senso).
  static std::vector<double> evaluate_exact_magnitude(
      const Grid<dim> &grid, const Function<dim> &exact_solution) {

    std::vector<double> exact_per_cell(grid.getArea(), 0.0);

    if constexpr (dim == 2) {
      for (int y = 0; y < static_cast<int>(grid.size.y); ++y) {
        for (int x = 0; x < static_cast<int>(grid.size.x); ++x) {
          types::Coordinate<dim> coord{x, y};
          std::size_t idx = grid.scalar_index(coord);

          types::ValuePoint<dim> phys_point{static_cast<double>(x),
                                            static_cast<double>(y)};

          utils::Vector<double, dim> u_exact = exact_solution.value(phys_point);
          exact_per_cell[idx] =
              std::sqrt(u_exact.dx * u_exact.dx + u_exact.dy * u_exact.dy);
        }
      }
    } else {
      static_assert(assertion::always_false<dim>,
                    "ErrorEvaluator::evaluate_exact_magnitude() non ancora implementato per dim > 2!");
    }

    return exact_per_cell;
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

/**
 * Errore globale rispetto a una soluzione analitica generica, con
 * norma scelta dal chiamante (default L2). Riusa ErrorEvaluator sia
 * per l'errore che per il termine di normalizzazione, così relative
 * = ||u_sim - u_exact|| / ||u_exact|| è calcolato con la stessa norma
 * su entrambi i lati.
 *
 * NOTA sulla semantica di L2_squared: il rapporto relative in quel
 * caso è tra somme di quadrati, non tra norme -- resta adimensionale
 * e monotono, ma non è confrontabile numericamente con un relative
 * calcolato in L2, L1 o Linfty.
 */
template <unsigned short int dim>
inline NormErrorResult compute_error(
    const Grid<dim> &grid,
    const std::vector<utils::Vector<double, dim>> &sim_u,
    const Function<dim> &exact_solution,
    NormType norm_type = NormType::L2) {

  const auto error_per_cell =
      ErrorEvaluator<dim>::integrate_difference(grid, sim_u, exact_solution);
  const auto exact_per_cell =
      ErrorEvaluator<dim>::evaluate_exact_magnitude(grid, exact_solution);

  const double abs_err =
      ErrorEvaluator<dim>::compute_global_error(error_per_cell, norm_type);
  const double ref_norm =
      ErrorEvaluator<dim>::compute_global_error(exact_per_cell, norm_type);

  if (ref_norm == 0.0)
    throw std::runtime_error(
        "compute_error(): reference norm is zero, cannot normalize");

  return NormErrorResult{abs_err / ref_norm, abs_err, norm_type};
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_ERROR_HPP