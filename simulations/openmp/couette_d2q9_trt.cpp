// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/simulation-config.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
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

  config::SimulationConfig<DIM> cfg;

  cfg.backend = lbm::ExecutionBackend::OPEN_MP;
  cfg.collision = lbm::TRT;

  cfg.grid_size = {129, 129};

  cfg.u0 = {0.1, 0};

  cfg.reynolds = 100;

  cfg.niters = 100000;
  cfg.nframes = 200;
  cfg.frames_out = "output/coutette_trt_frames";
  cfg.profile_out = "output/coutette_trt_profile.dat";
  // --- 2. ISTANZIA LOGGER ------------------------------------------------
  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  const DimPoint<DIM> grid_size(cfg.grid_size);

  LBM_LOG_INFO(
      main_logger,
      "Simulation:\n\tGrid dimensions: {}\n\tReynolds number: "
      "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
      "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
      grid_size, cfg.reynolds, cfg.u0, cfg.niters, cfg.nframes, cfg.frames_out,
      cfg.profile_out);

  // --- 3. CREA OSTACOLI --------------------------------------------------
  // Couette: parete inferiore rigida, parete superiore mobile, lati sinistro
  // e destro periodici. Gli angoli restano alle orizzontali come prima: il
  // wrap su x avviene per primo, poi la faccia y rivendica il link.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::PERIODIC;        // x = 0
  dbc.high(0) = Solid::PERIODIC;       // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;   // y = 0
  dbc.high(1) = Solid::BB_MOVING_WALL; // y = ny-1

  // --- 4. CREA MASCHERA --------------------------------------------------
  // Nessun ostacolo immerso nel fluido: la maschera e' tutta types::FLUID.
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out);

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc,
      CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, cfg.u0));
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_y_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  // H = altezza canale (parete inferiore a y=0, superiore a y=ny-1);
  // Umax = velocita' di riferimento (parete mobile per Couette).
  const double H = static_cast<double>(grid_size.y - 1);
  const auto exact_solution = analysis::CouetteSolution2D(H, cfg.u0.dx);
  const double err_l2 =
      simulation.compute_error(analysis::NormType::L2, exact_solution);

  LBM_LOG_NOTICE(main_logger, "{} error: {}",
                 analysis::to_string(analysis::NormType::L2), err_l2);

  simulation.detachListener(writer);
  solver.detachListener(writer);

  return 0;
}
