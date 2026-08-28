// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries/boundary-conditions.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm/config/config-parser.hpp"
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
#include <vector>

static constexpr lbm::types::dim_t DIM = 2;

/// Tipo di problema di cui si occupa questo binario.
static constexpr const char *PROBLEM_TYPE = "poiseuille";

/// L'operatore di collisione e' un parametro *template*, quindi va fissato a
/// compile time: questo binario esegue solo le configurazioni TRT, quelle
/// BGK le prende poiseuille_flow_2d_bgk.
static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::TRT;

/// Stessa logica per il backend.
static constexpr const char *BACKEND = "openmp";

namespace {

void print_usage(const char *exe) {
  std::cerr
      << "Uso: " << exe << " <config.toml> <output_vtk> <output_txt>\n\n"
      << "  config.toml   configurazione della simulazione\n"
      << "  output_vtk    CARTELLA in cui viene scritta la serie ParaView:\n"
      << "                <nome_config>_00000.vti, ... e il .pvd da aprire\n"
      << "                in ParaView. Viene creata se non esiste.\n"
      << "  output_txt    file col profilo di velocita' lungo la centerline\n";
}

} // namespace

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  // --- 1. LEGGI CONFIGURAZIONE ------------------------------------------
  if (argc != 4) {
    print_usage(argv[0]);
    return 1;
  }

  config::SimulationConfig cfg;
  try {
    cfg = config::parse_config(argv[1]);
    config::ensure_compatible(cfg, PROBLEM_TYPE, COLLISION, BACKEND, DIM);
  } catch (const config::ConfigError &err) {
    std::cerr << "Errore di configurazione: " << err.what() << "\n";
    return 1;
  }

  cfg.frames_out = argv[2];
  cfg.profile_out = argv[3];

  // --- 2. ISTANZIA LOGGER ------------------------------------------------
  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  const DimPoint<DIM> grid_size = cfg.grid_size<DIM>();
  const utils::Vector<double, DIM> init_vel = cfg.velocity<DIM>();

  LOG_INFO(main_logger,
           "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
           "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
           "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
           cfg.name, grid_size, cfg.reynolds, init_vel, cfg.niters,
           cfg.nframes, cfg.frames_out, cfg.profile_out);

  // --- 3. CREA OSTACOLI --------------------------------------------------
  // Poiseuille: pareti rigide sopra e sotto, ingresso e uscita a pressione
  // imposta sui lati. Gli angoli restano alle orizzontali come prima: il wrap
  // su x avviene per primo, poi la faccia y rivendica il link.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::PRESSURE_PERIODIC_INLET;   // x = 0
  dbc.high(0) = Solid::PRESSURE_PERIODIC_OUTLET; // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;             // y = 0
  dbc.high(1) = Solid::BB_RIGID_WALL;            // y = ny-1

  // --- 4. CREA MASCHERA --------------------------------------------------
  // Nessun ostacolo immerso nel fluido: la maschera e' tutta types::FLUID.
  types::solid_mask_t solid_mask = Solid::compute_solid_mask<DIM>({}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size,
                                               init_vel);
  // Salto di pressione che sostiene il flusso: ricavato dalla soluzione di
  // Poiseuille per il canale, non e' un parametro libero.
  const double pout = 1;
  const double pin =
      pout + 
      numbers::invcs_2 *
      (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) *
                 8 * params.nu * params.init_vel.dx;

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc, params, pin, pout);
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_y_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  const double H = static_cast<double>(grid_size.y - 1);
  const auto exact_solution = analysis::PoiseuilleSolution2D(H, init_vel.dx);
  const double err_l2 =
      simulation.compute_error(analysis::NormType::L2, exact_solution);

  LOG_NOTICE(main_logger, "{} error: {}",
             analysis::to_string(analysis::NormType::L2), err_l2);

  simulation.detachListener(writer);
  solver.detachListener(writer);

  return 0;
}
