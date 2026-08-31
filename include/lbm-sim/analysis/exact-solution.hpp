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
 * @section coords Coordinate convention
 *
 * A Coordinate<dim> passed to value() is a *node index*, not a physical
 * position. With halfway bounce-back the wall sits midway between the last
 * fluid node and the first solid node, so for a channel spanned by N fluid
 * rows the walls are at y = -1/2 and y = N - 1/2, and the physical distance
 * of node j from the bottom wall is
 *
 *     y_phys(j) = j + 1/2,      channel height H = N.
 *
 * Every solution below applies that +1/2 shift itself, via
 * @c halfway_wall_offset. Consequently the @c channel_height passed to the
 * constructors is the *number of fluid nodes across the channel*
 * (grid.size.y), never N - 1. Getting this half cell wrong is the usual
 * reason a converged run still reports a large error: it leaves an O(1/N)
 * discrepancy concentrated at the near-wall nodes, which on a 129-row grid
 * is about 1.5% of Umax at the first node -- large enough to swamp any real
 * difference between collision operators.
 */

#pragma once

#include "lbm-sim/functions.hpp"
#include "lbm-sim/logging.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <utility>
#include <vector>

namespace lbm {
namespace analysis {

/**
 * Distance, in lattice units, between a boundary node and the wall it
 * carries. Halfway bounce-back puts the wall half a spacing outside the last
 * fluid node; a full-way (on-node) bounce-back would use 0.0 instead.
 */
inline constexpr double halfway_wall_offset = 0.5;

/**
 * Analytical solution for 2D Couette flow.
 * Fixed bottom wall (y = -1/2), moving top wall (y = H - 1/2) with velocity
 * Umax. Linear velocity profile: u_x(y) = Umax * (y_phys / H), u_y = 0,
 * with y_phys = node index + halfway_wall_offset.
 *
 * @param channel_height  number of fluid rows across the channel (grid.size.y)
 * @param max_velocity    velocity of the moving wall
 */
class CouetteSolution2D : public functional::Function<2> {
private:
  double H;      // Channel height == number of fluid rows
  double Umax;   // Velocity of the moving top wall
  double offset; // Wall offset (halfway bounce-back by default)

public:
  static constexpr std::string_view name = "couette";

  CouetteSolution2D(double channel_height, double max_velocity,
                    double wall_offset = halfway_wall_offset)
      : H(channel_height), Umax(max_velocity), offset(wall_offset) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    // p.y is a node index; the wall lies `offset` below node 0.
    const double y_phys = static_cast<double>(p.y) + offset;
    const double ux = Umax * (y_phys / H);
    return utils::Vector<double, 2>{ux, 0.0};
  }

  const std::string_view &getName() const override { return name; }
};

/**
 * Analytical solution for 2D Poiseuille flow.
 * Driven by a pressure gradient / body force between two fixed walls at
 * y = -1/2 and y = H - 1/2. Parabolic velocity profile:
 * u_x(y) = 4 * Umax * (y_phys / H) * (1 - y_phys / H), u_y = 0,
 * with y_phys = node index + halfway_wall_offset.
 *
 * With this convention the first fluid node is not zero but
 * 4 * Umax * (0.5/H) * (1 - 0.5/H), i.e. about 0.0154 * Umax for H = 129.
 *
 * @param channel_height          number of fluid rows across the channel
 * @param max_center_velocity     peak velocity at the channel centre
 */
class PoiseuilleSolution2D : public functional::Function<2> {
private:
  double H;      // Channel height == number of fluid rows
  double Umax;   // Peak velocity at the channel centre
  double offset; // Wall offset (halfway bounce-back by default)

public:
  static constexpr std::string_view name = "poiseuille";

  PoiseuilleSolution2D(double channel_height, double max_center_velocity,
                       double wall_offset = halfway_wall_offset)
      : H(channel_height), Umax(max_center_velocity), offset(wall_offset) {}

  utils::Vector<double, 2> value(const types::Coordinate<2> &p) const override {
    const double y_phys = static_cast<double>(p.y) + offset;
    const double y_norm = y_phys / H;
    const double ux = 4.0 * Umax * y_norm * (1.0 - y_norm);
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
 * R is the *effective* wall radius, i.e. r_inner + halfway_wall_offset, the
 * same half-cell shift the plane cases apply along y. Prefer the
 * from_inner_radius() factory over passing R by hand.
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

  /**
   * Builds the solution from the radius handed to CylindricalShell, applying
   * the halfway bounce-back offset for you.
   */
  static HagenPoiseuilleSolution3D
  from_inner_radius(double r_inner, double max_axis_velocity, double axis_y,
                    double axis_z, double wall_offset = halfway_wall_offset) {
    return HagenPoiseuilleSolution3D(r_inner + wall_offset, max_axis_velocity,
                                     axis_y, axis_z);
  }

  utils::Vector<double, 3> value(const types::Coordinate<3> &p) const override {
    const double dy = static_cast<double>(p.y) - cy;
    const double dz = static_cast<double>(p.z) - cz;
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

/**
 * Writes a (coordinate, value) profile of an exact solution to disk.
 *
 * The coordinate column is the wall-normal position normalised by the
 * channel height, using the same halfway convention as the solutions:
 * (node index + halfway_wall_offset) / N. The value column is divided by
 * u_ref, so the simulated profile it is compared against must be normalised
 * the same way, or u_ref must be 1.0 here.
 */
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

/**
 * Normalised wall-normal position of a node, matching the convention the
 * solutions use. n_nodes is the number of fluid nodes across the channel.
 */
inline double normalised_wall_normal_coord(int index, int n_nodes) {
  return (static_cast<double>(index) + halfway_wall_offset) /
         static_cast<double>(n_nodes);
}

template <types::dim_t dim>
inline std::vector<std::pair<double, double>>
extract_dx_profile_along_y_center(const functional::Function<dim> &solution,
                                  const Grid<dim> &grid) {
  const int N = static_cast<int>(grid.size.y);
  const int x = grid.size.x / 2;
  std::vector<std::pair<double, double>> profile(N);
  for (int y = 0; y < N; ++y) {
    const double coord = normalised_wall_normal_coord(y, N);
    if constexpr (dim == 2) {
      profile[y] = {coord, solution.value({x, y}).dx};
    } else {
      profile[y] = {coord, solution.value({x, y, grid.size.z / 2}).dx};
    }
  }
  return profile;
}

template <types::dim_t dim>
inline std::vector<std::pair<double, double>>
extract_dy_profile_along_x_center(const functional::Function<dim> &solution,
                                  const Grid<dim> &grid) {
  const int N = static_cast<int>(grid.size.x);
  const int y = grid.size.y / 2;
  std::vector<std::pair<double, double>> profile(N);

  for (int x = 0; x < N; ++x) {
    const double coord = normalised_wall_normal_coord(x, N);
    if constexpr (dim == 2) {
      profile[x] = {coord, solution.value({x, y}).dy};
    } else {
      profile[x] = {coord, solution.value({x, y, grid.size.z / 2}).dy};
    }
  }

  return profile;
}

inline std::vector<std::pair<double, double>>
extract_dx_profile_along_z_center(const functional::Function<3> &solution,
                                  const Grid<3> &grid) {
  const int N = static_cast<int>(grid.size.z);
  const int x = grid.size.x / 2;
  const int y = grid.size.y / 2;
  std::vector<std::pair<double, double>> profile(N);

  for (int z = 0; z < N; ++z) {
    const double coord = normalised_wall_normal_coord(z, N);
    profile[z] = {coord, solution.value({x, y, z}).dx};
  }
  return profile;
}

} // namespace analysis
} // namespace lbm
