// Flusso in un tubo: Hagen-Poiseuille 3D in un condotto a sezione circolare,
// velocity set D3Q19, operatore BGK, backend OpenMP.
//
// E' il fratello 3D di poiseuille_flow_2d_bgk.cpp: stesso motore (pressione
// imposta su ingresso e uscita, bounce-back sulla parete), ma il canale
// piano fra due pareti diventa un condotto cilindrico. Cambia di
// conseguenza il profilo atteso, parabolico di rivoluzione invece che
// parabolico in y, e cambia il salto di pressione che lo sostiene:
//
//   canale (2D):  u_max = dp * H^2 / (8 * mu * L)
//   tubo   (3D):  u_max = dp * R^2 / (4 * mu * L)      <- Hagen-Poiseuille
//
// La parete non e' una faccia del dominio ma un ostacolo immerso: il
// dominio resta un parallelepipedo e la CylindricalShell dichiara solido
// tutto cio' che sta oltre il raggio del tubo, scavando il condotto dentro
// la scatola. Gli spigoli della scatola restano solidi e non partecipano.

// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/simulation-config.hpp"
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

static constexpr lbm::types::dim_t DIM = 3;
constexpr auto COLLISION = lbm::CollisionModel::BGK;
constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;

  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  // --- 1. CONFIGURAZIONE -------------------------------------------------
  config::SimulationConfig<DIM> cfg;

  cfg.name = "pipe_flow_3d_bgk";

  cfg.backend = BACKEND;
  cfg.collision = COLLISION;

  // Tubo lungo x. La sezione e' quadrata e il condotto ci sta dentro con
  // due celle di margine per lato: la parete cade cosi' *dentro* il
  // dominio, dove il bounce-back la puo' rappresentare, e non sulle facce.
  cfg.grid_size = {128, 65, 65};

  cfg.u0 = {0.1, 0, 0};

  cfg.reynolds = 100;

  cfg.niters = 100000;
  cfg.nframes = 200;

  cfg.frames_out = "output/pipe_flow_3d_bgk_frames";
  cfg.profile_out = "output/pipe_flow_3d_bgk_profile.dat";

  // --- 2. ISTANZIA LOGGER ------------------------------------------------
  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  const DimPoint<DIM> grid_size(cfg.grid_size);

  const int nx = static_cast<int>(grid_size.x);
  const int ny = static_cast<int>(grid_size.y);
  const int nz = static_cast<int>(grid_size.z);

  // Asse del tubo: la retta (y, z) = (cy, cz), parallela a x.
  const int cy = ny / 2;
  const int cz = nz / 2;
  const unsigned int radius = 30;

  LOG_INFO(main_logger,
           "Simulation '{}':\n\tGrid dimensions: {}\n\tPipe axis: (y,z) = "
           "({},{}), radius {}\n\tReynolds number: {}\n\tReference velocity: "
           "{}\n\tNumber of Iterations: {}\n\tNumber of frames: {}\n\tFrames "
           "output: {}\n\tProfile output: {}",
           cfg.name, grid_size, cy, cz, radius, cfg.reynolds, cfg.u0,
           cfg.niters, cfg.nframes, cfg.frames_out, cfg.profile_out);

  // --- 3. CREA OSTACOLI --------------------------------------------------
  // Ingresso e uscita a pressione imposta sulle due facce x; le altre
  // quattro facce sono gia' sepolte nel solido della parete, quindi la BC
  // che portano non viene mai raggiunta da un nodo di fluido: rigida e'
  // comunque la scelta coerente.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::PRESSURE_PERIODIC_INLET;   // x = 0
  dbc.high(0) = Solid::PRESSURE_PERIODIC_OUTLET; // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;             // y = 0
  dbc.high(1) = Solid::BB_RIGID_WALL;            // y = ny-1
  dbc.low(2) = Solid::BB_RIGID_WALL;             // z = 0
  dbc.high(2) = Solid::BB_RIGID_WALL;            // z = nz-1

  // La parete del tubo: guscio cilindrico da x = 0 a x = nx-1 attorno
  // all'asse. Il raggio esterno (ny + nz) non e' un raggio vero, e' solo
  // "abbastanza grande da uscire dal dominio": la AABB viene tosata sulla
  // griglia, quindi il guscio copre di fatto ogni nodo oltre `radius`,
  // spigoli della scatola compresi.
  const CollisionDetection::CollisionArea<DIM> pipe_wall(
      Coordinate<DIM>(0, 0, 0),
      {CollisionDetection::CylindricalShell<DIM>(Coordinate<DIM>(0, cy, cz),
                                                 Coordinate<DIM>(nx - 1, cy, cz),
                                                 radius, ny + nz)});

  // id 0 = la parete del tubo: bounce-back rigido, ferma.
  const std::vector<Solid::ObstacleData<DIM>> obstacle_data{
      {Solid::BB_RIGID_WALL, {0.0, 0.0, 0.0}}};

  // --- 4. CREA MASCHERA --------------------------------------------------
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({pipe_wall}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  // ATTENZIONE: CollisionParams ricava nu dalla lunghezza caratteristica
  // num_cells.y, non dal diametro del tubo. Con ny = 65 e un diametro
  // effettivo 2R+1 = 61 il Reynolds davvero simulato e' u*D/nu ~ 94, non
  // 100: tenere il raggio vicino a ny/2 se si vuole che i due coincidano.
  const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size, cfg.u0);

  // Salto di pressione che sostiene il flusso: come nel canale 2D non e' un
  // parametro libero, lo fissa la soluzione di Hagen-Poiseuille invertita
  // per dare u_max = init_vel.dx sull'asse. In unita' di reticolo p = rho*cs^2,
  // quindi il salto di *densita'* e' invcs_2 volte quello di pressione.
  //
  // R e' il raggio effettivo della parete: col bounce-back halfway la parete
  // sta a meta' strada fra l'ultimo nodo di fluido e il primo nodo solido.
  const double R_eff = radius + 0.5;
  const double pout = 1;
  const double pin = pout + numbers::invcs_2 * 4.0 * params.nu *
                                params.init_vel.dx * nx / (R_eff * R_eff);

  LBMSimulation<DIM, D3Q19, COLLISION> simulation(
      grid_size, std::move(solid_mask), obstacle_data, dbc, params, pin, pout);
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D3Q19, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  // ux lungo z sulla colonna centrale (x = nx/2, y = ny/2): e' un diametro
  // del tubo, cioe' la sezione su cui si legge la parabola.
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_z_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  const auto exact_solution =
      analysis::HagenPoiseuilleSolution3D(R_eff, cfg.u0.dx, cy, cz);
  const double err_l2 =
      simulation.compute_error(analysis::NormType::L2, exact_solution);

  LOG_NOTICE(main_logger, "{} error: {}",
             analysis::to_string(analysis::NormType::L2), err_l2);

  simulation.detachListener(writer);
  solver.detachListener(writer);

  return 0;
}
