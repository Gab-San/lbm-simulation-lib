#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/formatting.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"
#include "lbm/config/config-parser.hpp"

// Path to the Ghia benchmarks, injected by CMake (see
// simulations/CMakeLists.txt) so it does not depend on the working directory.
// The fallback only applies when building outside CMake.
#ifndef LBM_BENCHMARKS_DIR
#define LBM_BENCHMARKS_DIR "benchmarks"
#endif

// C++ STD LIB
static constexpr unsigned short int DIM = 2;

/// Lid-driven cavity: three rigid walls + the moving lid on top.
static lbm::Solid::DomainBC<DIM> make_cavity_bc() {
  lbm::Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = lbm::Solid::BB_RIGID_WALL;   // x = 0
  dbc.high(0) = lbm::Solid::BB_RIGID_WALL;  // x = nx-1
  dbc.low(1) = lbm::Solid::BB_RIGID_WALL;   // y = 0
  dbc.high(1) = lbm::Solid::BB_MOVING_WALL; // the lid, y = ny-1
  return dbc;
}

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;

  if (argc < 2) {
    config::print_usage(argv[0]);
    return 1;
  }

  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  std::vector<config::SimulationConfig<DIM>> configs;
  try {
    configs = config::parse_config<DIM>(argv[1]);
  } catch (const config::ConfigError &err) {
    LBM_LOG_CRITICAL(main_logger, "Config Error {}", err.what());
    return 1;
  }

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  const std::string path_to_benchmark =
      std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);
    const utils::Vector<double, DIM> u0(cfg.u0);

    LBM_LOG_INFO(
        main_logger,
        "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
        "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
        "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
        cfg.name, grid_size, cfg.reynolds, u0, cfg.niters, cfg.nframes,
        cfg.frames_out, cfg.profile_out);

    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>({}, grid_size);

    std::shared_ptr<VtkWriter> writer =
        std::make_shared<VtkWriter>(cfg.frames_out);

    Simulation simulation(
        grid_size, std::move(solid_mask), {}, make_cavity_bc(),
        CollisionParams<DIM, CollisionType>(cfg.reynolds, grid_size, u0));
    simulation.attachListener(writer);

    CUDASolver<DIM, D2Q9, CollisionType> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/);

    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dy_profile_along_x_center<DIM>);

    // Comparison against Ghia et al. (1982). Norm chosen here: L2.
    const auto ghia_y = simulation.compute_ghia_error(
        path_to_benchmark + "data_y_" + format::file_format(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_y.relative,
                   ghia_y.absolute);

    const auto ghia_x = simulation.compute_ghia_error(
        path_to_benchmark + "data_x_" + format::file_format(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_x.relative,
                   ghia_x.absolute);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
