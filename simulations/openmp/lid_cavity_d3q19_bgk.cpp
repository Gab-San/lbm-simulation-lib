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

static constexpr lbm::types::dim_t DIM = 3;

static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

namespace {

/// Maschera dei confini per una cavita' Nx*Ny*Nz: 5 pareti rigide
/// (bounce-back semplice) + 1 parete mobile (il "lid", bounce-back con
/// velocita' imposta) sulla faccia superiore z = Nz-1.
///
/// Non serve piu' scorrere le facce nodo per nodo: sono esattamente le sei
/// facce del dominio, sei byte in tutto.
lbm::Solid::DomainBC<3> build_domain_bc() {
  lbm::Solid::DomainBC<3> dbc{};
  for (lbm::types::dim_t a = 0; a < 3; ++a) {
    dbc.low(a) = lbm::Solid::BB_RIGID_WALL;
    dbc.high(a) = lbm::Solid::BB_RIGID_WALL;
  }
  dbc.high(2) = lbm::Solid::BB_MOVING_WALL; // il lid, z = nz-1
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
    std::cerr << "Errore di configurazione: " << err.what() << "\n";
    return 1;
  }
  // --- 2. ISTANZIA LOGGER ------------------------------------------------
  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");
  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);

  const utils::Vector<double, DIM> lid_velocity = cfg.u0;

  // ATTENZIONE: in 3D il numero di nodi cresce come N^3 e il costo per nodo
  // e' 19/9 volte quello di D2Q9. Una cavita' 128^3 ha gia' oltre 2M di
  // nodi: conviene validare su 32^3 o 48^3 prima di scalare.
  LBM_LOG_INFO(main_logger,
               "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
               "{}\n\tLid velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
               "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
               cfg.name, grid_size, cfg.reynolds, lid_velocity, cfg.niters,
               cfg.nframes, cfg.frames_out, cfg.profile_out);

  // --- 3./4. CREA OSTACOLI E MASCHERA ------------------------------------
  // Le pareti sono le facce del dominio; nessun ostacolo immerso nel fluido,
  // quindi la maschera e' tutta types::FLUID.
  const Solid::DomainBC<3> dbc = build_domain_bc();
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  //
  // Lo stesso writer va agganciato sia a `sim` sia a `solver`: il primo
  // notifica l'header con le dimensioni della griglia, il secondo i frame
  // delle norme di velocita', e sono due DataObservable distinti.
std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

  LBMSimulation<DIM, D3Q19, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc,
      CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, lid_velocity));
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D3Q19, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_z_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  // Non c'e': la lid cavity 3D non ha soluzione analitica, e le tabelle di
  // Ghia et al. sono solo per il caso 2D (compute_ghia_error() e' infatti
  // definita solo per dim == 2). Il confronto quantitativo va fatto contro
  // dati di riferimento esterni sul profilo scritto sopra.

  simulation.detachListener(writer);
  solver.detachListener(writer);

   #ifdef LBM_PROFILING
  lbm::profiling::dump_csv(cfg.profile_out);  
  lbm::profiling::reset();                    
#endif
  }
  return 0;
}
