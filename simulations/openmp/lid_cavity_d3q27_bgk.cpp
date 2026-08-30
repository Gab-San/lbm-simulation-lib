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
#include <iostream>
#include <memory>
#include <string>

static constexpr lbm::types::dim_t DIM = 3;

static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

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
    const utils::Vector<double, DIM> lid_velocity = cfg.u0;

    // WARNING: in 3D the node count grows as N^3 and the per-node cost is
    // 19/9 times that of D2Q9. A 128^3 cavity already has more than 2M nodes:
    // validate on 32^3 or 48^3 before scaling up.
    LBM_LOG_INFO(main_logger,
                 "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
                 "{}\n\tLid velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
                 "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
                 cfg.name, grid_size, cfg.reynolds, lid_velocity, cfg.niters,
                 cfg.nframes, cfg.frames_out, cfg.profile_out);

    // --- 3./4. CREATE OBSTACLES AND MASK -----------------------------------
    // The walls are the domain faces; no obstacle is immersed in the fluid, so
    // the mask is entirely types::FLUID.
    const Solid::DomainBC<3> dbc = build_domain_bc();
    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>({}, grid_size);

    // --- 5. RUN SIMULATION -------------------------------------------------
    // frames_out is the DIRECTORY; the file basename comes from the config
    // name, so different runs in the same directory do not overwrite each
    // other.
    //
    // The same writer has to be attached both to `sim` and to `solver`: the
    // former notifies the header with the grid dimensions, the latter the
    // frames of velocity norms, and they are two distinct DataObservable.
    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

    LBMSimulation<DIM, D3Q27, COLLISION> simulation(
        grid_size, std::move(solid_mask), {}, dbc,
        CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, lid_velocity));
    simulation.attachListener(writer);

    OpenMPSolver<DIM, D3Q27, COLLISION> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver);

    // --- 6. OUTPUT ---------------------------------------------------------
    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dx_profile_along_z_center);

    // --- 7. ERROR COMPUTATION ----------------------------------------------
    // There is none: the 3D lid cavity has no analytical solution, and the
    // Ghia et al. tables only cover the 2D case (compute_ghia_error() is in
    // fact defined only for dim == 2). The quantitative comparison has to be
    // made against external reference data, on the profile written above.

    simulation.detachListener(writer);
    solver.detachListener(writer);

#ifdef LBM_PROFILING
    lbm::profiling::dump_csv(cfg.profile_out);
    lbm::profiling::reset();
#endif
  }
  return 0;
}
