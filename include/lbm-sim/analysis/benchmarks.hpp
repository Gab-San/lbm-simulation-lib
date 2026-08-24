#ifndef __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP
#define __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP

// LBM SIM LIB
#include "lbm-sim/analysis/error.hpp"
#include "lbm-sim/analysis/types.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/lattice.hpp"

// C++ STANDARD LIB
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace lbm {
namespace formatting {

inline std::string format_reyn(double reyn) {
  std::ostringstream oss;
  oss << reyn;
  return oss.str();
}

} // namespace formatting

namespace analysis {

namespace detail {

struct GhiaCavityData {
  static constexpr std::size_t n_points = 17;
  std::array<double, n_points> coord;
  std::array<double, n_points> values;

  GhiaCavityData(const std::array<double, n_points> &coord_,
                 const std::array<double, n_points> &values_)
      : coord(coord_), values(values_) {};
  ~GhiaCavityData() = default;
};

inline void check_path_existence(const std::string &filepath_in) {
  std::filesystem::path path_in(filepath_in);
  if (!path_in.has_filename()) {
    throw std::runtime_error("GHIA_INPUT: Filepath does not point to a file!");
  }

  if (!std::filesystem::exists(path_in)) {
    throw std::runtime_error("GHIA_INPUT: File not found!");
  }
}

inline std::vector<std::string> split(std::string s,
                                      const std::string &delimiter) {
  std::vector<std::string> tokens;
  size_t pos = 0;

  std::string token;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    tokens.push_back(token);
    s.erase(0, pos + delimiter.length());
  }

  tokens.push_back(s);
  return tokens;
}

enum PARSING_TYPE { Y, X };

// Assumes that the benchmark file is correctly formatted.
inline enum PARSING_TYPE parse_header(std::ifstream &in, double reyn) {
  std::string header;

  // FINISH READING HEADER
  std::getline(in, header);
  const std::vector<std::string> header_toks = split(header, " ");

  if (header_toks[0].substr(0, 2) != "%%") {
    throw std::runtime_error(
        "(Badly Formatted Header) Error while parsing GHIA_INPUT: missing %%");
  }

  if (reyn != std::stod(header_toks[3])) {
    std::string error =
        "Error while parsing GHIA_INPUT: mismatching reynolds number " +
        std::to_string(reyn) + "|" + std::to_string(std::stod(header_toks[3]));

    throw std::runtime_error(error);
  }

  return header_toks[2] == "y" ? Y : X;
}

inline GhiaCavityData read_data(std::ifstream &input) {
  std::array<double, GhiaCavityData::n_points> coord;
  std::array<double, GhiaCavityData::n_points> values;

  auto lineidx = 0;
  for (std::string line; std::getline(input, line);) {
    std::istringstream iss(line);
    double line_coord, line_value;

    // This automatically skips all whitespace and extracts the doubles
    if (iss >> line_coord >> line_value) {
      coord[lineidx] = line_coord;
      values[lineidx] = line_value;
      lineidx++;
    }
  }
  if (lineidx != GhiaCavityData::n_points) {
    throw std::runtime_error("Error while parsing GHIA_INPUT: number of rows "
                             "is less than expected!");
  }

  return GhiaCavityData(coord, values);
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

inline NormErrorResult
compute_ghia_error(const std::string &filepath_in, const Lattice<2> &lattice,
                   const double reynolds, const double lid_velocity,
                   const NormType norm_type = NormType::L2) {
  using namespace detail;
  check_path_existence(filepath_in);

  std::ifstream in(filepath_in);
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open ghia benchmark file!");
  }

  enum PARSING_TYPE pt = detail::parse_header(in, reynolds);

  GhiaCavityData data = detail::read_data(in);

  // Extract simulated u along vertical or horizontal centerline,
  // normalized by lid_velocity, so it's directly comparable to Ghia's
  // tabulated values.
  const bool is_x = (pt == PARSING_TYPE::X);

  const auto line_size = is_x ? lattice.grid.size.x : lattice.grid.size.y;
  const auto fixed_coord =
      (is_x ? lattice.grid.size.y : lattice.grid.size.x) / 2;

  std::vector<double> velocity_along_coord_line(line_size);

  for (unsigned int i = 0; i < line_size; ++i) {
    const types::Coordinate<2> p = is_x ? types::Coordinate<2>(i, fixed_coord)
                                        : types::Coordinate<2>(fixed_coord, i);

    const auto &vel = lattice.u[lattice.grid.scalar_index(p)];
    const double component = is_x ? vel.dy : vel.dx;

    velocity_along_coord_line[i] = component / lid_velocity;
  }

  // error_per_point / ref_per_point contain magnitudes (>= 0):
  // ErrorEvaluator::compute_global_error() squares them or takes the maximum
  // value based on norm_type, therefore we can pass
  // |diff| and |ref value|, not the signed values.
  std::array<double, detail::GhiaCavityData::n_points> u_diff_per_point;
  std::array<double, detail::GhiaCavityData::n_points> u_ref_per_point;

  for (auto i = 0; i < detail::GhiaCavityData::n_points; ++i) {
    const double u_sim =
        detail::interp_along_line(velocity_along_coord_line, data.coord[i]);
    u_diff_per_point[i] = std::abs(u_sim - data.values[i]);
    u_ref_per_point[i] = std::abs(data.values[i]);
  }

  const double u_abs_err =
      ErrorEvaluator<2>::compute_global_error(u_diff_per_point, norm_type);
  const double u_ref_norm =
      ErrorEvaluator<2>::compute_global_error(u_ref_per_point, norm_type);

  if (u_ref_norm == 0.0) {
    throw std::runtime_error(
        "compute_ghia_error(): reference norm is zero, check the table");
  }

  return NormErrorResult{u_abs_err / u_ref_norm, u_abs_err, norm_type};
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_GHIA_BENCHMARK_HPP
