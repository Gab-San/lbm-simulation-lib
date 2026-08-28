// Variante di lid_cavity_d2q9_bgk che gira il solve dentro uno
// scopedApply() di BackendProperties: serve a misurare l'effetto di un
// numero di thread imposto senza toccare il resto della pipeline.

#include "lbm-sim/backend/properties.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/config/config-parser.hpp"
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
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
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

/// Numero di thread imposto per la durata del solve: e' il parametro sotto
/// test di questo binario.
static constexpr int NUM_THREADS = 6;

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;

  config::SimulationConfig<DIM> cfg;

  cfg.name = "properties_test";

  cfg.backend = lbm::ExecutionBackend::OPEN_MP;
  cfg.collision = lbm::CollisionModel::BGK;

  cfg.grid_size = {129, 129};

  cfg.u0 = {0.1, 0};

  cfg.reynolds = 100;

  cfg.niters = 100000;
  cfg.nframes = 200;

  cfg.frames_out = "output/properties_test_frames";
  cfg.profile_out = "output/properties_test_profile.dat";

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

  auto &prop = profiling::BackendProperties<OPEN_MP>::get();
  prop.setNumThreads(NUM_THREADS);

  // --- 3. CREA OSTACOLI --------------------------------------------------
  // Per la lid cavity gli "ostacoli" sono i quattro lati del dominio: tre
  // pareti rigide + il lid mobile in alto. Ricavati da grid_size.
  const int x_max = static_cast<int>(cfg.grid_size[0]) - 1;
  const int y_max = static_cast<int>(cfg.grid_size[1]) - 1;

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, y_max);
  const Coordinate<2> C(x_max, y_max);
  const Coordinate<2> D(x_max, 0);

  const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles{
      CollisionDetection::CollisionArea(A, {CollisionDetection::Segment(A, B),
                                            CollisionDetection::Segment(A, D),
                                            CollisionDetection::Segment(D, C)}),
      CollisionDetection::CollisionArea(A,
                                        {CollisionDetection::Segment(B, C)})};

  const std::unordered_map<unsigned int, uint8_t> obst_type_map{
      {0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}};

  // --- 4. CREA MASCHERA --------------------------------------------------
  types::boundary_mask_t obstacle_mask =
      Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, obstacle_mask,
      CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, init_vel));
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  {
    // Le proprieta' valgono solo per la durata del solve: alla chiusura
    // dello scope il runtime OpenMP torna com'era, cosi' una run non
    // condiziona quella dopo.
    const auto scope = prop.scopedApply();
    simulation.solve(solver);
  }

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
