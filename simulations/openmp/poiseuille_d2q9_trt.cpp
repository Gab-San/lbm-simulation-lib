// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/config-parser.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"

// C++ STD LIB
#include <memory>
#include <string>
#include <vector>

static constexpr lbm::types::dim_t DIM = 2;
constexpr auto COLLISION = lbm::CollisionModel::TRT;
constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;

  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

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

  LBM_LOG_INFO(
      main_logger,
      "Simulation:\n\tGrid dimensions: {}\n\tReynolds number: "
      "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
      "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
      grid_size, cfg.reynolds, cfg.u0, cfg.niters, cfg.nframes, cfg.frames_out,
      cfg.profile_out);

  // --- 3. CREATE OBSTACLES -----------------------------------------------
  // Poiseuille: rigid top and bottom walls, pressure-imposed inlet and outlet
  // on the sides. Corners still belong to the horizontal faces as before: the
  // wrap on x happens first, then the y face claims the link.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::PRESSURE_PERIODIC_INLET;   // x = 0
  dbc.high(0) = Solid::PRESSURE_PERIODIC_OUTLET; // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;             // y = 0
  dbc.high(1) = Solid::BB_RIGID_WALL;            // y = ny-1

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

  const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size, cfg.u0);
  // The pressure drop that sustains the flow: derived from the channel
  // Poiseuille solution, it is not a free parameter.
  const double pout = 1;
  const double pin =
      pout +
      numbers::invcs_2 *
          (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) * 8 *
          params.nu * params.init_vel.dx;

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc, params, pin, pout);
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_y_center);

  // --- 7. ERROR COMPUTATION ----------------------------------------------
  // H = channel height (bottom wall at y=0, top wall at y=ny-1);
  // Umax = reference velocity (the moving wall, for Couette).
  const double H = static_cast<double>(grid_size.y - 1);
  const auto exact_solution = analysis::PoiseuilleSolution2D(H, cfg.u0.dx);
  const double err_l2 =
      simulation.compute_error(analysis::NormType::L2, exact_solution);

  LBM_LOG_NOTICE(main_logger, "{} error: {}",
                 analysis::to_string(analysis::NormType::L2), err_l2);

  simulation.detachListener(writer);
  solver.detachListener(writer);
           #ifdef LBM_PROFILING
  lbm::profiling::dump_csv(cfg.profile_out);  
  lbm::profiling::reset();                    
#endif
}
  return 0;
}
