/**
 * @file exact-solution.hpp
 * @brief The closed-form solutions the library validates against.
 *
 * Three canonical cases, all steady-state and all in lattice units, so they
 * can be compared with @c lattice.u directly and without rescaling:
 * plane Couette flow, plane Poiseuille flow, and Hagen-Poiseuille flow in a
 * cylindrical pipe.
 *
 * Each is a Function<dim>, so it is passed straight to
 * LBMSimulation::compute_error(). A case with no closed form -- the 3D lid
 * cavity -- has no entry here and is validated against external reference
 * data instead.
 *
 * @note The wall position assumed by these solutions is the halfway
 *       bounce-back one: the wall sits midway between the last fluid node
 *       and the first solid node, not on the solid nodes. Getting that half
 *       cell wrong is the usual reason a converged run still reports a large
 *       error.
 */

#ifndef __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP
#define __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP

#include "lbm-sim/functions.hpp"
#include "lbm-sim/logging.hpp"

#include <filesystem>
#include <fstream>
#include <functional>

namespace lbm {
namespace analysis {

/**
 * Analytical solution for 2D Couette flow.
 * Fixed bottom wall (y = 0), moving top wall (y = H) with velocity Umax.
 * Linear velocity profile: u_x(y) = Umax * (y / H), u_y = 0.
 */
class CouetteSolution2D : public functional::Function<2> {
private:
  double H;    // Channel height
  double Umax; // Velocity of the moving top wall

public:
  static constexpr std::string_view name = "couette";
  CouetteSolution2D(double channel_height, double max_velocity)
      : H(channel_height), Umax(max_velocity) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    // p.y is the vertical coordinate.
    double ux = Umax * (p.y / H);
    return utils::Vector<double, 2>{ux, 0.0};
  }

  const std::string_view &getName() const override { return name; }
};

/**
 * Analytical solution for 2D Poiseuille flow.
 * Driven by a pressure gradient / body force between two fixed walls
 * (y = 0 and y = H). Parabolic velocity profile:
 * u_x(y) = 4 * Umax * (y / H) * (1 - y / H), u_y = 0.
 */
class PoiseuilleSolution2D : public functional::Function<2> {
private:
  double H;    // Channel height
  double Umax; // Peak velocity at the channel centre (y = H/2)

public:
  static constexpr std::string_view name = "poiseuille";

  PoiseuilleSolution2D(double channel_height, double max_center_velocity)
      : H(channel_height), Umax(max_center_velocity) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    double y_norm = p.y / H;
    double ux = 4.0 * Umax * y_norm * (1.0 - y_norm);
    return utils::Vector<double, 2>{ux, 0.0};
  }

  const std::string_view &getName() const override { return name; }
};

/**
 * Analytical solution for Hagen-Poiseuille flow in a 3D cylindrical pipe
 * whose axis is parallel to x.
 * Parabolic profile of revolution: u_x(r) = Umax * (1 - r^2/R^2), with r the
 * distance from the pipe axis; u_y = u_z = 0.
 *
 * Outside the pipe (r >= R) it evaluates to zero rather than to the negative
 * continuation of the parabola: that way the solid wall nodes, where the
 * solver leaves u = 0, do not pollute the error that
 * ErrorEvaluator<3>::integrate_difference() computes over the whole grid.
 *
 * R is the *effective* wall radius. With halfway bounce-back the wall does
 * not sit on the solid nodes but midway between the last fluid node and the
 * first solid node, hence R = r_inner + 0.5 where r_inner is the radius
 * passed to CylindricalShell.
 */
class HagenPoiseuilleSolution3D : public functional::Function<3> {
private:
  double R;      // Effective pipe radius
  double Umax;   // Velocity on the axis
  double cy, cz; // Axis position within the cross-section plane

public:
  static constexpr std::string_view name = "hagen-poiseuille";

  HagenPoiseuilleSolution3D(double pipe_radius, double max_axis_velocity,
                            double axis_y, double axis_z)
      : R(pipe_radius), Umax(max_axis_velocity), cy(axis_y), cz(axis_z) {}

  utils::Vector<double, 3> value(const types::Coordinate<3> &p) const override {
    const double dy = p.y - cy;
    const double dz = p.z - cz;
    const double r2 = dy * dy + dz * dz;
    const double R2 = R * R;

    if (r2 >= R2) {
      return utils::Vector<double, 3>{0.0, 0.0, 0.0};
    }

    const double ux = Umax * (1.0 - r2 / R2);
    return utils::Vector<double, 3>{ux, 0.0, 0.0};
  }

  const std::string_view &getName() const override { return name; }
};

template <types::dim_t dim, typename ExtractFn>
void dump_exact_solution_points(const std::filesystem::path &filepath,
                                const Grid<dim> &grid,
                                const functional::Function<dim> &exact_solution,
                                ExtractFn &&extract_profile, double u_ref) {
  logging::Logger *data_logger = logging::create_or_get_logger("data_log");

  std::ofstream fout(filepath, std::ios::binary);
  if (!fout.is_open()) {
    LBM_LOG_ERROR(
        data_logger,
        "Failed to create file: {}\nExact solution will not be dumped!",
        filepath.string());
    return;
  }

  const std::vector<std::pair<double, double>> profile =
      extract_profile(exact_solution, grid);

  std::vector<double> flat;
  flat.reserve(profile.size() * 2);
  for (const auto &[coord, value] : profile) {
    flat.push_back(coord);
    flat.push_back(value / u_ref);
  }

  const std::string header = "%%exact " +
                             std::string(exact_solution.getName()) + " " +
                             std::to_string(profile.size()) + "\n";

  fout.write(header.data(), static_cast<std::streamsize>(header.size()));
  fout.write(reinterpret_cast<const char *>(flat.data()),
             static_cast<std::streamsize>(flat.size() * sizeof(double)));

  if (!fout) {
    LBM_LOG_ERROR(data_logger, "Write failed for file: {}", filepath.string());
    return;
  }
  LBM_LOG_DEBUG(data_logger, "[File: {}] Profile generation complete...",
                filepath.string());
}

template <types::dim_t dim>
inline std::vector<std::pair<double, double>>
extract_dx_profile_along_y_center(const functional::Function<dim> &solution,
                                  const Grid<dim> &grid) {
  const int N = static_cast<int>(grid.size.y);
  const int x = grid.size.x / 2;
  std::vector<std::pair<double, double>> profile(N);
  for (int y = 0; y < N; ++y) {
    const double coord = (y + 0.5) / N;
    if constexpr (dim == 2) {
      profile[y] = {coord, solution.value({x, y}).dx};
    } else {
      profile[y] = {coord, solution.value({x, y, grid.size.z / 2}).dx};
    }
  }
  return profile;
}

template <types::dim_t dim>
inline std::vector<std::pair<float, float>>
extract_dy_profile_along_x_center(const functional::Function<dim> &solution,
                                  const Grid<dim> &grid) {
  const int N = static_cast<int>(grid.size.y);
  std::vector<std::pair<float, float>> profile(N);
  const int y = grid.size.y / 2;

  for (int x = 0; x < N; ++x) {
    const double coord = (x + 0.5) / N;
    if constexpr (dim == 2) {
      profile[x] = {coord, solution.value({x, y}).dy};
    } else {
      profile[x] = {coord, solution.value({x, y, grid.size.z / 2}).dy};
    }
  }

  return profile;
}

inline std::vector<std::pair<float, float>>
extract_dx_profile_along_z_center(const functional::Function<3> &solution,
                                  const Grid<3> &grid) {
  const int N = static_cast<int>(grid.size.z);
  std::vector<std::pair<float, float>> profile(N);
  const int x = grid.size.x / 2;
  const int y = grid.size.y / 2;

  for (int z = 0; z < N; ++z) {
    const double coord = (z + 0.5) / N;
    profile[z] = {coord, solution.value({x, y, z}).dx};
  }
  return profile;
}

} // namespace analysis
} // namespace lbm

#endif // __LBM_SIM_ANALYSIS_EXACT_SOLUTIONS_HPP
