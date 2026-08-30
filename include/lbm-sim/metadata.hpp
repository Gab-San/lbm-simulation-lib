/**
 * @file metadata.hpp
 * @brief Compile-time tags that parameterise the library, plus their
 *        human-readable names.
 *
 * Both enumerations are unscoped, so they convert implicitly to @c int and
 * can be used directly as non-type template arguments (@c cm_t of
 * LBMSimulation, @c backend_t of SolverBase). Selecting the collision model
 * and the backend at compile time keeps virtual dispatch out of the inner
 * loop.
 */

#pragma once

#include <string>

namespace lbm {

/**
 * @brief Execution path of a solver.
 *
 * Used as the @c backend_t template parameter of SolverBase and deduced by
 * LBMSimulation::solve() from the solver it is handed.
 */
enum ExecutionBackend {
  CUDA,   ///< GPU path (optional, requires the CUDA build).
  OPEN_MP ///< Multi-threaded CPU path.
};

/**
 * @brief Collision operator used to relax the distributions.
 *
 * Also selects the matching CollisionParams<dim, cm_t> specialisation, so
 * the parameter set is checked at compile time against the model.
 */
enum CollisionModel {
  BGK, ///< Single-relaxation-time (Bhatnagar-Gross-Krook).
  TRT, ///< Two-relaxation-time: symmetric/antisymmetric splitting.
  MRT  ///< Multiple-relaxation-time: each moment relaxed separately.
};

/**
 * @brief Display name of a collision model.
 *
 * @param cm_t Collision model to convert.
 * @return "BGK", "TRT" or "MRT".
 *
 * @warning The trailing @c std::to_string() is a fallback for values
 *          outside the enumeration and is unreachable for the enumerators
 *          above. A new enumerator added without its @c case will silently
 *          stringify as a number instead of failing to compile.
 *
 * @see LBMSimulation::output(), which stamps this name into the profile
 *      file header.
 */
inline std::string collision_model_to_string(enum CollisionModel cm_t) {
  switch (cm_t) {
  case CollisionModel::BGK:
    return "BGK";
  case CollisionModel::TRT:
    return "TRT";
  case CollisionModel::MRT:
    return "MRT";
  }
  return std::to_string(cm_t);
}

/**
 * @brief Display name of an execution backend.
 *
 * @param backend_t Backend to convert.
 * @return "OpenMP" or "CUDA".
 *
 * @warning Same fallback caveat as collision_model_to_string().
 */
inline std::string backend_to_string(enum ExecutionBackend backend_t) {
  switch (backend_t) {
  case ExecutionBackend::OPEN_MP:
    return "OpenMP";
  case ExecutionBackend::CUDA:
    return "CUDA";
  }
  return std::to_string(backend_t);
}

} // namespace lbm
