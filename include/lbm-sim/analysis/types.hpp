/**
 * @file types.hpp
 * @brief The vocabulary of the analysis layer: an abstract exact solution,
 *        the norms, and the result of a comparison.
 *
 * Deliberately thin and free of any dependency on the solver. Anything that
 * can answer "what is the exact velocity at this point" is a Function<dim>,
 * so a user-defined analytical solution plugs into
 * LBMSimulation::compute_error() with no change to the library.
 */

#pragma once

#include <cstddef>

namespace lbm::analysis {

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
 * Outcome of a comparison against a generic analytical solution
 * (Couette, Poiseuille, ...), with the norm chosen by the caller.
 */
struct NormErrorResult {
  double relative;
  double absolute;
  NormType norm_type;

  double rmse = 0.0;
  double linf = 0.0;
  double rmse_normalized = 0.0;
  double linf_normalized = 0.0;
  std::size_t sample_count = 0;
};

} // namespace lbm::analysis
