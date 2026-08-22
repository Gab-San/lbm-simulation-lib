#ifndef __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP
#define __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP

// LBM SIM LIB
#include "lbm-sim/analysis/error.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/lattice.hpp"

// C++ STANDARD LIB
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lbm {
namespace analysis {

struct GhiaCavityData {
  static constexpr std::size_t n_points = 17;

  // y-coordinate (normalized, 0 = bottom wall, 1 = moving lid) for
  // Table I, u-velocity along the vertical centerline x = 0.5.
  static constexpr std::array<double, n_points> y = {
      0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531,
      0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766,
      1.0000};

  static constexpr std::array<double, n_points> u_re100 = {
      0.00000,  -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
      -0.15662, -0.21090, -0.20581, -0.13641, 0.00332,  0.23151,
      0.68717,  0.73722,  0.78871,  0.84123,  1.00000};

  static constexpr std::array<double, n_points> u_re400 = {
      0.00000,  -0.08186, -0.09266, -0.10338, -0.14612, -0.24299,
      -0.32726, -0.17119, -0.11477, 0.02135,  0.16256,  0.29093,
      0.55892,  0.61756,  0.68439,  0.75837,  1.00000};

  static constexpr std::array<double, n_points> u_re1000 = {
      0.00000,  -0.18109, -0.20196, -0.22220, -0.29730, -0.38289,
      -0.27805, -0.10648, -0.06080, 0.05702,  0.18719,  0.33304,
      0.46604,  0.51117,  0.57492,  0.65928,  1.00000};

  // x-coordinate (normalized, 0 = left wall, 1 = right wall) for
  // Table II, v-velocity along the horizontal centerline y = 0.5.
  static constexpr std::array<double, n_points> x = {
      0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344,
      0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688,
      1.0000};

  static constexpr std::array<double, n_points> v_re100 = {
      0.00000, 0.09233,  0.10091,  0.10890,  0.12317,  0.16077,
      0.17507, 0.17527,  0.05454,  -0.24533, -0.22445, -0.16914,
      -0.10313, -0.08864, -0.07391, -0.05906, 0.00000};

  // NOTE: entry at x=0.9063 (-0.23827) flagged as possibly wrong by
  // the transcriber -- verify against the original paper before
  // relying on it.
  static constexpr std::array<double, n_points> v_re400 = {
      0.00000, 0.18360,  0.19713,  0.20920,  0.22965,  0.28124,
      0.30203, 0.30174,  0.05186,  -0.38598, -0.44993, -0.23827,
      -0.22847, -0.19254, -0.15663, -0.12146, 0.00000};

  static constexpr std::array<double, n_points> v_re1000 = {
      0.00000, 0.27485,  0.29012,  0.30353,  0.32627,  0.37095,
      0.33075, 0.32235,  0.02526,  -0.31966, -0.42665, -0.51500,
      -0.39188, -0.33714, -0.27669, -0.21388, 0.00000};
};

enum class GhiaReynolds { Re100, Re400, Re1000 };

/**
 * Converte un numero di Reynolds "grezzo" (es. da Config::reyn_num) nel
 * GhiaReynolds corrispondente. Lancia se il valore non e' uno dei tre
 * per cui abbiamo dati tabulati (100, 400, 1000) -- niente arrotondamenti
 * silenziosi verso il Reynolds piu' vicino.
 */
inline GhiaReynolds ghia_reynolds_from_value(double reyn_num) {
  if (reyn_num == 100.0)
    return GhiaReynolds::Re100;
  if (reyn_num == 400.0)
    return GhiaReynolds::Re400;
  if (reyn_num == 1000.0)
    return GhiaReynolds::Re1000;
  throw std::invalid_argument(
      "ghia_reynolds_from_value(): no Ghia reference data for this Reynolds "
      "number (only 100, 400, 1000 are tabulated)");
}

namespace detail {

inline const std::array<double, GhiaCavityData::n_points> &
u_table_for(GhiaReynolds re) {
  switch (re) {
  case GhiaReynolds::Re100:
    return GhiaCavityData::u_re100;
  case GhiaReynolds::Re400:
    return GhiaCavityData::u_re400;
  case GhiaReynolds::Re1000:
    return GhiaCavityData::u_re1000;
  }
  throw std::invalid_argument("u_table_for(): unhandled GhiaReynolds value");
}

inline const std::array<double, GhiaCavityData::n_points> &
v_table_for(GhiaReynolds re) {
  switch (re) {
  case GhiaReynolds::Re100:
    return GhiaCavityData::v_re100;
  case GhiaReynolds::Re400:
    return GhiaCavityData::v_re400;
  case GhiaReynolds::Re1000:
    return GhiaCavityData::v_re1000;
  }
  throw std::invalid_argument("v_table_for(): unhandled GhiaReynolds value");
}

// Linear interpolation of the simulated field at a normalized
// coordinate t in [0,1] along a grid line of length n_cells.
inline double interp_along_line(const std::vector<double> &line_values,
                                 double t) {
  const double pos = t * static_cast<double>(line_values.size() - 1);
  const std::size_t i0 = static_cast<std::size_t>(std::floor(pos));
  const std::size_t i1 = std::min(i0 + 1, line_values.size() - 1);
  const double frac = pos - static_cast<double>(i0);
  return line_values[i0] * (1.0 - frac) + line_values[i1] * frac;
}

} // namespace detail

/**
 * Compares the simulated u-velocity along the vertical centerline
 * (x = Nx/2) against Ghia's Table I, and the simulated v-velocity
 * along the horizontal centerline (y = Ny/2) against Ghia's Table II.
 *
 * lid_velocity: the top-wall speed used to drive the cavity (Ghia's
 * values are normalized by it -- typically params.init_vel.dx).
 *
 * norm_type: quale norma usare per ridurre i 17 punti tabulati a un
 * singolo scalare (default L2), stessa semantica di NormType usata in
 * error.hpp / compute_error(). L'errore assoluto risultante e' sempre
 * calcolato con la STESSA norma usata per il termine di normalizzazione
 * (relative = ||diff||_norm / ||ref||_norm), cosi' come per compute_error.
 *
 * Returns one NormErrorResult for u and one for v, each computed over
 * the 17 tabulated Ghia points (simulated field linearly interpolated
 * to those normalized positions, not the other way around).
 */
struct GhiaComparisonResult {
  NormErrorResult u_error;
  NormErrorResult v_error;
};

inline GhiaComparisonResult compute_ghia_error(const Lattice<2> &lattice,
                                                GhiaReynolds re,
                                                double lid_velocity,
                                                NormType norm_type = NormType::L2) {
  if (lid_velocity == 0.0) {
    throw std::invalid_argument(
        "compute_ghia_error(): lid_velocity must be nonzero");
  }

  const auto &u_ref = detail::u_table_for(re);
  const auto &v_ref = detail::v_table_for(re);

  // Extract simulated u along vertical centerline (x fixed, y varies),
  // normalized by lid_velocity, so it's directly comparable to Ghia's
  // tabulated values.
  const unsigned int x_c = lattice.grid.size.x / 2;
  std::vector<double> u_line(lattice.grid.size.y);
  for (unsigned int y = 0; y < lattice.grid.size.y; ++y) {
    const types::Coordinate<2> p(x_c, y);
    u_line[y] = lattice.u[lattice.grid.scalar_index(p)].dx / lid_velocity;
  }

  // Extract simulated v along horizontal centerline (y fixed, x varies).
  const unsigned int y_c = lattice.grid.size.y / 2;
  std::vector<double> v_line(lattice.grid.size.x);
  for (unsigned int x = 0; x < lattice.grid.size.x; ++x) {
    const types::Coordinate<2> p(x, y_c);
    v_line[x] = lattice.u[lattice.grid.scalar_index(p)].dy / lid_velocity;
  }

  // error_per_point / ref_per_point contengono magnitudini (>= 0):
  // ErrorEvaluator::compute_global_error() eleva al quadrato/prende il
  // massimo internamente a seconda di norm_type, quindi qui passiamo
  // |diff| e |valore di riferimento|, non i valori con segno.
  std::vector<double> u_diff_per_point(GhiaCavityData::n_points);
  std::vector<double> u_ref_per_point(GhiaCavityData::n_points);
  for (std::size_t i = 0; i < GhiaCavityData::n_points; ++i) {
    const double u_sim = detail::interp_along_line(u_line, GhiaCavityData::y[i]);
    u_diff_per_point[i] = std::abs(u_sim - u_ref[i]);
    u_ref_per_point[i] = std::abs(u_ref[i]);
  }

  std::vector<double> v_diff_per_point(GhiaCavityData::n_points);
  std::vector<double> v_ref_per_point(GhiaCavityData::n_points);
  for (std::size_t i = 0; i < GhiaCavityData::n_points; ++i) {
    const double v_sim = detail::interp_along_line(v_line, GhiaCavityData::x[i]);
    v_diff_per_point[i] = std::abs(v_sim - v_ref[i]);
    v_ref_per_point[i] = std::abs(v_ref[i]);
  }

  const double u_abs_err =
      ErrorEvaluator<2>::compute_global_error(u_diff_per_point, norm_type);
  const double u_ref_norm =
      ErrorEvaluator<2>::compute_global_error(u_ref_per_point, norm_type);

  const double v_abs_err =
      ErrorEvaluator<2>::compute_global_error(v_diff_per_point, norm_type);
  const double v_ref_norm =
      ErrorEvaluator<2>::compute_global_error(v_ref_per_point, norm_type);

  if (u_ref_norm == 0.0 || v_ref_norm == 0.0) {
    throw std::runtime_error(
        "compute_ghia_error(): reference norm is zero, check the table");
  }

  return GhiaComparisonResult{
      NormErrorResult{u_abs_err / u_ref_norm, u_abs_err, norm_type},
      NormErrorResult{v_abs_err / v_ref_norm, v_abs_err, norm_type}};
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP