#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/config-parser.hpp"
#include "lbm-sim/config/simulation-config.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"
#include "lbm/logging.hpp"

// QUILL LIB
#include "quill/LogMacros.h"

// C++ STD LIB
#include <memory>
#include <string>
#include <vector>

// Path ai benchmark di Ghia, iniettato da CMake (vedi
// simulations/CMakeLists.txt) per non dipendere dalla working directory.
// Il fallback vale solo se si compila fuori da CMake.
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

  config::SimulationConfig<DIM> cfg;

  cfg.name = "lid_cavity_d2q9_bgk";

  cfg.backend = lbm::ExecutionBackend::OPEN_MP;
  cfg.collision = lbm::CollisionModel::BGK;

  cfg.grid_size = {129, 129};

  cfg.u0 = {0.1, 0};

  cfg.reynolds = 100;

  cfg.niters = 100000;
  cfg.nframes = 200;

  cfg.frames_out = "output/lid_cavity_bgk_frames";
  cfg.profile_out = "output/lid_cavity_bgk_profile.dat";

  // --- 2. ISTANZIA LOGGER ------------------------------------------------
  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  const DimPoint<DIM> grid_size(cfg.grid_size);
  const utils::Vector<double, DIM> init_vel = cfg.u0;

  LOG_INFO(main_logger,
           "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
           "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
           "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
           cfg.name, grid_size, cfg.reynolds, init_vel, cfg.niters, cfg.nframes,
           cfg.frames_out, cfg.profile_out);

  // --- 3. CREA OSTACOLI --------------------------------------------------

  // Per la lid cavity i "confini" non sono ostacoli disegnati sui nodi di
  // bordo, ma le quattro facce del dominio: tre pareti rigide + il lid mobile
  // in alto. Quattro byte in tutto, indipendenti dalla risoluzione.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::BB_RIGID_WALL;   // x = 0
  dbc.high(0) = Solid::BB_RIGID_WALL;  // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;   // y = 0
  dbc.high(1) = Solid::BB_MOVING_WALL; // il lid, y = ny-1

  // --- 4. CREA MASCHERA --------------------------------------------------
  // Nessun ostacolo immerso nel fluido: la maschera e' tutta types::FLUID.
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc,
      CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, init_vel));
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dy_profile_along_x_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  // Confronto con Ghia et al. (1982). Norma scelta qui: L2.
  const std::string path_to_benchmark =
      std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

  const auto ghia_y = simulation.compute_ghia_error(
      path_to_benchmark + "data_y_" + formatting::format_reyn(cfg.reynolds) +
      ".txt");

  LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
             analysis::to_string(analysis::NormType::L2), ghia_y.relative,
             ghia_y.absolute);

  const auto ghia_x = simulation.compute_ghia_error(
      path_to_benchmark + "data_x_" + formatting::format_reyn(cfg.reynolds) +
      ".txt");

  LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
             analysis::to_string(analysis::NormType::L2), ghia_x.relative,
             ghia_x.absolute);

  simulation.detachListener(writer);
  solver.detachListener(writer);

  return 0;
}
