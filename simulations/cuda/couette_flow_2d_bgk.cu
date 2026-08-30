// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries/boundary-conditions.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/config-parser.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"

// C++ STD LIB
#include <memory>
#include <vector>

static constexpr unsigned short int DIM = 2;

/// Couette: rigid bottom wall, moving top wall, periodic sides.
static lbm::Solid::DomainBC<DIM> make_couette_bc() {
  lbm::Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = lbm::Solid::PERIODIC;        // x = 0
  dbc.high(0) = lbm::Solid::PERIODIC;       // x = nx-1
  dbc.low(1) = lbm::Solid::BB_RIGID_WALL;   // y = 0
  dbc.high(1) = lbm::Solid::BB_MOVING_WALL; // y = ny-1
  return dbc;
}

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

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

  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);

    LBM_LOG_INFO(
        main_logger,
        "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
        "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
        "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
        cfg.name, grid_size, cfg.reynolds, cfg.u0, cfg.niters, cfg.nframes,
        cfg.frames_out, cfg.profile_out);

    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>({}, grid_size);

    std::shared_ptr<VtkWriter> writer =
        std::make_shared<VtkWriter>(cfg.frames_out);

    Simulation simulation(
        grid_size, std::move(solid_mask), {}, make_couette_bc(),
        CollisionParams<DIM, CollisionType>(cfg.reynolds, grid_size, cfg.u0));
    simulation.attachListener(writer);

    CUDASolver<DIM, D2Q9, CollisionType> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver);
    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dx_profile_along_y_center);

    // H = channel height (bottom wall at y=0, top at y=grid_size.y-1);
    // Umax = reference velocity (the moving wall, for Couette).
    // The same values already used to build the simulation: no duplication,
    // flow_type picks the right Function<2>.
    const double H = static_cast<double>(grid_size.y - 1);
    const auto exact_solution = analysis::CouetteSolution2D(H, cfg.u0.dx);
    const double err_l2 =
        simulation.compute_error(analysis::NormType::L2, exact_solution);

    LBM_LOG_NOTICE(main_logger, "{} error: {}",
                   analysis::to_string(analysis::NormType::L2), err_l2);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
