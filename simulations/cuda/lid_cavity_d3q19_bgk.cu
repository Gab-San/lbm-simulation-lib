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

// C++ STD LIB
static constexpr unsigned short int DIM = 3;

namespace {

/// Boundary mask for an Nx*Ny*Nz cavity: 5 rigid walls (plain bounce-back)
/// + 1 moving wall (the "lid", bounce-back with imposed velocity) on the top
/// face z = Nz-1.
///
/// There is no need to walk the faces node by node any more: they are exactly
/// the six domain faces, six bytes in total.
lbm::Solid::DomainBC<3> build_domain_bc() {
  lbm::Solid::DomainBC<3> dbc{};
  for (lbm::types::dim_t a = 0; a < 3; ++a) {
    dbc.low(a) = lbm::Solid::BB_RIGID_WALL;
    dbc.high(a) = lbm::Solid::BB_RIGID_WALL;
  }
  dbc.high(2) = lbm::Solid::BB_MOVING_WALL; // the lid, z = nz-1
  return dbc;
}

} // namespace

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
  using Simulation = LBMSimulation<DIM, D3Q19, CollisionType>;

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
        grid_size, std::move(solid_mask), {}, build_domain_bc(),
        CollisionParams<DIM, CollisionType>(cfg.reynolds, grid_size, u0));
    simulation.attachListener(writer);

    CUDASolver<DIM, D3Q19, CollisionType> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);
    simulation.solve(solver /*, preconditioner*/);

    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dx_profile_along_z_center);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
