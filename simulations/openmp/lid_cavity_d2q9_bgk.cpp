// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/formatting.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"
#include "lbm/config/config-parser.hpp"

// C++ STD LIB
#include <memory>
#include <string>
#include <vector>

// Path to the Ghia benchmarks, injected by CMake (see
// simulations/CMakeLists.txt) so it does not depend on the working directory.
// The fallback only applies when building outside CMake.
#ifndef LBM_BENCHMARKS_DIR
#define LBM_BENCHMARKS_DIR "benchmarks"
#endif

static constexpr lbm::types::dim_t DIM = 2;

static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;

  if (argc < 2) {
    config::print_usage(argv[0]);
    return 1;
  }

  std::vector<config::SimulationConfig<DIM>> configs;
  try {
    configs = config::parse_config<DIM>(argv[1]);
  } catch (const config::ConfigError &err) {
    std::cerr << "Configuration error: " << err.what() << "\n";
    return 1;
  }
  // --- 2. INSTANTIATE LOGGER ---------------------------------------------
  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);

    const utils::Vector<double, DIM> init_vel = cfg.u0;

    LBM_LOG_INFO(
        main_logger,
        "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
        "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
        "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
        cfg.name, grid_size, cfg.reynolds, init_vel, cfg.niters, cfg.nframes,
        cfg.frames_out, cfg.profile_out);

    // --- 3. CREATE OBSTACLES -----------------------------------------------

    // For the lid cavity the "boundaries" are not obstacles painted on the
    // border nodes, but the four domain faces: three rigid walls + the moving
    // lid on top. Four bytes in total, independent of the resolution.
    Solid::DomainBC<DIM> dbc{};
    dbc.low(0) = Solid::BB_RIGID_WALL;   // x = 0
    dbc.high(0) = Solid::BB_RIGID_WALL;  // x = nx-1
    dbc.low(1) = Solid::BB_RIGID_WALL;   // y = 0
    dbc.high(1) = Solid::BB_MOVING_WALL; // the lid, y = ny-1

    // --- 4. CREATE MASK ----------------------------------------------------
    // No obstacle immersed in the fluid: the mask is entirely types::FLUID.
    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>({}, grid_size);

    // --- 5. RUN SIMULATION -------------------------------------------------
    // frames_out is the DIRECTORY; the file basename comes from the config
    // name, so different runs in the same directory do not overwrite each
    // other.
    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

    LBMSimulation<DIM, D2Q9, COLLISION> simulation(
        grid_size, std::move(solid_mask), {}, dbc,
        CollisionParams<DIM, COLLISION>(cfg.reynolds, init_vel, init_vel.dx,
                                        grid_size.y));
    simulation.attachListener(writer);

    OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver);

    // --- 6. OUTPUT ---------------------------------------------------------
    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dy_profile_along_x_center<DIM>);

    // --- 7. ERROR COMPUTATION ----------------------------------------------
    // Comparison against Ghia et al. (1982). Norm chosen here: L2.
    const std::string path_to_benchmark =
        std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

    const auto ghia_y = simulation.compute_ghia_error(
        path_to_benchmark + "data_y_" + format::file_format(cfg.reynolds) +
        ".txt");

    const auto ghia_x = simulation.compute_ghia_error(
        path_to_benchmark + "data_x_" + format::file_format(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger,
                   "Lid-driven cavity validation summary\n"
                   "  Grid:                 {} x {}\n"
                   "  Reynolds:             {}\n"
                   "  Collision:            {}\n"
                   "  Iterations:           {}\n\n"
                   "  Ghia ux(y/2):\n"
                   "    relative L2:        {:.4f} %\n"
                   "    RMSE / U_lid:       {:.4f} %\n"
                   "    Linf / U_lid:       {:.4f} %\n\n"
                   "  Ghia uy(x/2):\n"
                   "    relative L2:        {:.4f} %\n"
                   "    RMSE / U_lid:       {:.4f} %\n"
                   "    Linf / U_lid:       {:.4f} %\n\n",
                   grid_size.x, grid_size.y, cfg.reynolds,
                   collision_model_to_string(COLLISION), cfg.niters,
                   100.0 * ghia_x.relative, 100.0 * ghia_x.rmse_normalized,
                   100.0 * ghia_x.linf_normalized, 100.0 * ghia_y.relative,
                   100.0 * ghia_y.rmse_normalized,
                   100.0 * ghia_y.linf_normalized);

    format::CsvWriter<analysis::DetailedErrorAnalysisSchema> error_writer(
        "out/error_" + cfg.name + "_d2q9_" + format::file_format(grid_size) +
            "_" + format::file_format(cfg.reynolds) + "_" +
            format::file_format(COLLISION) + ".csv",
        true);

    error_writer.append_row("uy", grid_size.x, grid_size.y, cfg.reynolds,
                            collision_model_to_string(COLLISION), cfg.niters,
                            analysis::to_string(analysis::NormType::L2),
                            ghia_y.relative, ghia_y.absolute,
                            100.0 * ghia_y.relative, ghia_y.rmse,
                            100.0 * ghia_y.rmse_normalized, ghia_y.linf,
                            100.0 * ghia_y.linf_normalized);

    error_writer.append_row("ux", grid_size.x, grid_size.y, cfg.reynolds,
                            collision_model_to_string(COLLISION), cfg.niters,
                            analysis::to_string(analysis::NormType::L2),
                            ghia_x.relative, ghia_x.absolute,
                            100.0 * ghia_x.relative, ghia_x.rmse,
                            100.0 * ghia_x.rmse_normalized, ghia_x.linf,
                            100.0 * ghia_x.linf_normalized);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }
  return 0;
}
