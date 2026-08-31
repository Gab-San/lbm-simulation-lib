/**
 * @file lbm-simulation.hpp
 * @brief LBMSimulation: user-facing driver that owns the lattice and the
 *        collision parameters and runs a solver over them.
 */

#ifndef __LBM_SIM_LBM_SIMULATION_HPP
#define __LBM_SIM_LBM_SIMULATION_HPP

#include "lbm-sim/analysis/benchmarks.hpp"
#include "lbm-sim/analysis/error.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lattice.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/solver/solver-base.hpp"
#include "lbm-sim/types/common.hpp"

// C++ STANDARD LIB
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

namespace lbm {

/**
 * @brief Owns a Lattice and its CollisionParams, runs a solver on them and
 *        exposes the post-processing entry points.
 *
 * The caller never builds a Lattice directly: it passes the lattice's
 * constructor arguments through and the simulation constructs it in place.
 * After solve() has returned, the macroscopic fields inside the lattice are
 * what compute_error(), compute_ghia_error() and output() read.
 *
 * Deriving from DataObservable lets the simulation push raw byte buffers to
 * registered listeners; the grid header written by solve() is one such
 * buffer.
 *
 * @tparam dim         Spatial dimension (2 or 3).
 * @tparam VelocitySet Discrete velocity set (D2Q9, D3Q19, ...). Must expose
 *                     a @c dim member matching @p dim.
 * @tparam cm_t        Collision model; defaults to CollisionModel::BGK.
 */
template <unsigned short int dim, typename VelocitySet,
          enum CollisionModel cm_t = CollisionModel::BGK>
class LBMSimulation : public DataObservable {
private:
  static_assert(
      dim == VelocitySet::dim,
      "LBMSimulation: template parameter 'dim' must match VelocitySet::dim");

  /// Simulation state, owned by this object.
  Lattice<dim> lattice;

  /// Relaxation parameters, Reynolds number and initial velocity.
  const CollisionParams<dim, cm_t> params;

public:
  /**
   * @brief Builds the lattice in place and validates the domain boundary
   *        conditions.
   *
   * The first four arguments and the two pressures are forwarded to the
   * Lattice constructor; @p params_ is stored as-is. Consistency of the
   * domain boundary conditions (for instance a periodic side paired with a
   * non-periodic opposite side) is checked by
   * Solid::assert_consistent_domain_bc().
   *
   * @param grid_dim_   Domain extents.
   * @param solid_mask_ Per-cell obstacle ids (moved in).
   * @param obstacles_  Obstacle table indexed by id (moved in, may be empty).
   * @param domain_bc_  Boundary condition of each domain side.
   * @param params_     Collision parameters matching @p cm_t.
   * @param pin         Inlet pressure. Defaults to 0.
   * @param pout        Outlet pressure. Defaults to 0.
   */
  LBMSimulation(const types::DimPoint<dim> grid_dim_,
                types::solid_mask_t solid_mask_,
                std::vector<Solid::ObstacleData<dim>> obstacles_,
                Solid::DomainBC<dim> domain_bc_,
                const CollisionParams<dim, cm_t> params_, const double pin = 0,
                const double pout = 0)
      : lattice(grid_dim_, std::move(solid_mask_), std::move(obstacles_),
                domain_bc_, pin, pout),
        params(params_) {
    Solid::assert_consistent_domain_bc<dim>(domain_bc_);
  };

  /**
   * @brief Runs the simulation with the given solver.
   *
   * Logs the start, broadcasts the grid dimensions to the listeners via
   * write_header() so they know the shape of the data that follows, then
   * hands lattice and params to the solver. The solver evolves its own
   * distributions and writes the macroscopic fields back into
   * @c lattice.u and @c lattice.rho.
   *
   * @tparam backend_t Execution backend, deduced from @p solver, so the
   *                   same simulation type accepts either an OpenMP or a
   *                   CUDA solver.
   * @param solver     Solver to run. Taken by reference and not stored, so
   *                   the simulation may be solved again with a different
   *                   one -- bearing in mind the lattice still holds the
   *                   previous run's fields.
   */
  template <enum ExecutionBackend backend_t>
  void solve(SolverBase<dim, VelocitySet, cm_t, backend_t> &solver) {
    logging::Logger *simulation_logger =
        logging::create_or_get_logger("simulation");

    LBM_LOG_DEBUG(simulation_logger, "Initializing Simulation...");

    write_header(lattice.grid);

    solver.solve(lattice, params);

    LBM_LOG_DEBUG(simulation_logger, "Finished Simulation.");
  };

  /**
   * @brief Error against an analytical solution, in the style of
   *        dealii::VectorTools.
   *
   * The caller supplies only the exact solution (a Function<dim> such as
   * CouetteSolution2D / PoiseuilleSolution2D, or any user-defined derived
   * class) and the norm type; the approximate field (@c lattice.u) and the
   * grid are taken internally from the simulation and must not be passed
   * by the caller.
   *
   * Equivalent to:
   * @code
   *   VectorTools::integrate_difference(..., solution, exact_solution,
   *                                     error_per_cell, ..., norm_type);
   *   VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
   * @endcode
   *
   * @param norm_type      Norm used for the global reduction.
   * @param exact_solution Analytical solution to compare against.
   * @return The absolute (non-normalised) global error in the requested
   *         norm. For the error relative to the exact solution, see
   *         analysis::compute_error() in analysis/error.hpp.
   */
  double compute_error(const analysis::NormType &norm_type,
                       const functional::Function<dim> &exact_solution) const {
    const auto error_per_cell =
        analysis::ErrorEvaluator<dim>::integrate_difference(
            lattice.grid, lattice.u, exact_solution);

    const double error = analysis::ErrorEvaluator<dim>::compute_global_error(
        error_per_cell, norm_type);

    return error;
  }

  /**
   * @brief Error against the Ghia et al. (1982) benchmark, 2D lid-driven
   *        cavity only.
   *
   * Compares @c lattice.u along the two centrelines with the reference
   * tables. The lid velocity -- the velocity of the moving wall,
   * @c params.init_vel.dx -- is the one Ghia uses to normalise his data,
   * and @c params.reyn_num selects the tabulated Reynolds number.
   *
   * @param filepath_in Path to the Ghia reference tables.
   * @param norm_type   Norm used to reduce the 17 tabulated points to a
   *                    single scalar (default L2), same semantics as
   *                    compute_error().
   * @return The error in the requested norm.
   *
   * @note Available only for @c dim == 2: Ghia's lid cavity has no
   *       tabulated 3D equivalent in this library. Instantiating this
   *       member for @c dim == 3 is a compile error; a 3D simulation that
   *       never calls it compiles fine.
   */
  analysis::NormErrorResult compute_ghia_error(
      const std::string &filepath_in,
      const analysis::NormType norm_type = analysis::NormType::L2) const {
    if constexpr (dim == 2) {
      return analysis::compute_ghia_error(filepath_in, lattice, params.reyn_num,
                                          params.init_vel.dx);
    } else {
      static_assert(assertion::always_false<dim>,
                    "compute_ghia_error() is only defined for dim == 2");
    }
  }

  /**
   * @brief Extracts a 1D profile from the lattice and writes it to disk.
   *
   * Creates the parent directory chain if missing, then writes an ASCII
   * header line followed by the raw payload:
   * @verbatim
   * %%profile <MODEL> <N> <LID_VELOCITY>\n
   * <N doubles, native byte order>
   * @endverbatim
   * where @c MODEL is collision_model_to_string(cm_t), @c N is the profile
   * size and @c LID_VELOCITY is @c params.init_vel.dx.
   *
   * @param filepath        Destination path; parent directories are created.
   * @param extract_profile Callable reducing the lattice to a vector of
   *                        doubles, e.g. the velocity magnitude along a
   *                        centreline.
   *
   * @warning If the file cannot be opened the failure is logged and the
   *          function returns silently, so the caller cannot detect it
   *          (see the FIXME below).
   * @warning The payload is written verbatim, so the file is endianness-
   *          and double-representation-dependent: read it back on a
   *          matching platform.
   */
  void output(const char *filepath,
              std::function<std::vector<double>(const Lattice<dim> &)>
                  extract_profile) {
    using namespace std::filesystem;

    logging::Logger *data_logger = logging::create_or_get_logger("data_log");

    path outpath(filepath);
    path parent = outpath.parent_path();

    if (!exists(parent)) {
      create_directories(parent);
    }

    std::ofstream fout(outpath, std::ios::binary);

    if (!fout.is_open()) {
      LBM_LOG_ERROR(data_logger,
                    "Failed to create file: {}\nProfile for this simulation "
                    "will not be dumped!",
                    filepath);
      return;
    }

    std::vector<double> profile = extract_profile(lattice);

    LBM_LOG_DEBUG(data_logger, "Writing header to file {}", filepath);

    std::string header = "%%profile " + collision_model_to_string(cm_t) + " " +
                         std::to_string(profile.size()) + " " +
                         std::to_string(params.init_vel.dx) + "\n";

    fout.write(header.data(), header.size());

    LBM_LOG_DEBUG(data_logger, "Writing data to file {}", filepath);

    fout.write(reinterpret_cast<const char *>(profile.data()),
               profile.size() * sizeof(double));

    fout.close();

    LBM_LOG_DEBUG(data_logger, "[File: {}] Profile generation complete...",
                  filepath);
  };

private:
  /**
   * @brief Broadcasts the grid extents to the registered listeners.
   *
   * Packs @c (nx, ny) for @c dim == 2, @c (nx, ny, nz) otherwise, each cast
   * to @c int32_t so the wire format does not depend on the platform's
   * @c size_t width, and hands the buffer to notifyListeners(). Byte order
   * is native.
   *
   * Called once from solve() before the solver runs, so every listener
   * receives the domain shape before any field data.
   *
   * @param grid Grid whose extents are published.
   */
  void write_header(const Grid<dim> &grid) {
    std::vector<char> buf(sizeof(int32_t) * dim);

    if constexpr (dim == 2) {
      const int32_t nx32 = static_cast<int32_t>(grid.size.x);
      const int32_t ny32 = static_cast<int32_t>(grid.size.y);

      std::memcpy(buf.data(), &nx32, sizeof(int32_t));
      std::memcpy(buf.data() + sizeof(int32_t), &ny32, sizeof(int32_t));
    } else {
      const int32_t nx32 = static_cast<int32_t>(grid.size.x);
      const int32_t ny32 = static_cast<int32_t>(grid.size.y);
      const int32_t nz32 = static_cast<int32_t>(grid.size.z);

      std::memcpy(buf.data(), &nx32, sizeof(int32_t));
      std::memcpy(buf.data() + sizeof(int32_t), &ny32, sizeof(int32_t));
      std::memcpy(buf.data() + 2 * sizeof(int32_t), &nz32, sizeof(int32_t));
    }

    LBM_LOG_DEBUG(logging::create_or_get_logger("data_log"),
                  "Writing norms header...");

    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm

#endif // __LBM_SIM_LBM_SIMULATION_HPP
