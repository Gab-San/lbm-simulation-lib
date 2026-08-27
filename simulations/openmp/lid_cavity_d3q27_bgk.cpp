// Lid-driven cavity 3D, velocity set D3Q27, backend OpenMP.
//
// D3Q27 e' il set completo: 27 direzioni, incluse quelle di vertice.
// Piu' isotropo di D3Q19 ma piu' costoso per nodo (27 vs 19
// popolazioni da streammare e collidere).
//
// LBMSimulation e OpenMPSolver sono templati sul VelocitySet, quindi
// l'unica differenza col main D3Q19 e' il tipo passato ai due template.

#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm/config/config-parser.hpp"
#include "lbm-sim/core/grid.hpp"
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

static constexpr lbm::types::dim_t DIM = 3;

/// Tipo di problema di cui si occupa questo binario.
static constexpr const char *PROBLEM_TYPE = "lid_cavity";

/// L'operatore di collisione e' un parametro *template*, quindi va fissato a
/// compile time.
static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

/// Stessa logica per il backend.
static constexpr const char *BACKEND = "openmp";

namespace {

void print_usage(const char *exe) {
  std::cerr
      << "Uso: " << exe << " <config.toml> <output_vtk> <output_txt>\n\n"
      << "  config.toml   configurazione della simulazione. Deve essere 3D,\n"
      << "                cioe' dichiarare [grid].nz\n"
      << "  output_vtk    CARTELLA in cui viene scritta la serie ParaView:\n"
      << "                <nome_config>_00000.vti, ... e il .pvd da aprire\n"
      << "                in ParaView. Viene creata se non esiste.\n"
      << "  output_txt    file col profilo ux(z) lungo la verticale centrale\n";
}

/// Maschera dei confini per una cavita' Nx*Ny*Nz: 5 pareti rigide
/// (bounce-back semplice) + 1 parete mobile (il "lid", bounce-back con
/// velocita' imposta) sulla faccia superiore z = Nz-1.
///
/// A differenza del 2D non passiamo da Solid::compute_boundary_mask(): il
/// modulo collision-detection lavora su Segment/Circle nel piano (x,y) e
/// non e' ancora esteso a 3D. Per un dominio a scatola non serve comunque
/// quella generalita' - basta scorrere le facce.
lbm::types::boundary_mask_t
build_boundary_mask(const lbm::types::DimPoint<3> &size) {
  lbm::types::boundary_mask_t mask(size.x * size.y * size.z, lbm::Solid::NONE);
  const lbm::Grid<3> grid(size);

#pragma omp parallel for collapse(3) schedule(static)
  for (std::size_t z = 0; z < size.z; ++z) {
    for (std::size_t y = 0; y < size.y; ++y) {
      for (std::size_t x = 0; x < size.x; ++x) {
        const bool on_boundary =
            (x == 0 || x == size.x - 1 || y == 0 || y == size.y - 1 ||
             z == 0 || z == size.z - 1);
        if (!on_boundary)
          continue;

        const lbm::types::Coordinate<3> p(
            static_cast<int>(x), static_cast<int>(y), static_cast<int>(z));

        mask[grid.scalar_index(p)] = (z == size.z - 1)
                                         ? lbm::Solid::BB_MOVING_WALL
                                         : lbm::Solid::BB_RIGID_WALL;
      }
    }
  }

  return mask;
}

} // namespace

int main(int argc, char **argv) {
  using namespace lbm;
  using types::DimPoint;

  // --- 1. LEGGI CONFIGURAZIONE ------------------------------------------
  if (argc != 4) {
    print_usage(argv[0]);
    return 1;
  }

  config::SimulationConfig cfg;
  try {
    cfg = config::parse_config(argv[1]);
    // ensure_compatible controlla anche la dimensione: una config senza
    // [grid].nz viene letta come 2D e qui rifiutata, invece di girare
    // ignorando in silenzio la terza dimensione.
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
  const utils::Vector<double, DIM> lid_velocity = cfg.velocity<DIM>();

  // ATTENZIONE: in 3D il numero di nodi cresce come N^3 e il costo per nodo
  // e' 27/9 = 3x quello di D2Q9. Una cavita' 128^3 ha gia' oltre 2M di
  // nodi: conviene validare su 32^3 o 48^3 prima di scalare.
  LOG_INFO(main_logger,
           "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
           "{}\n\tLid velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
           "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
           cfg.name, grid_size, cfg.reynolds, lid_velocity, cfg.niters,
           cfg.nframes, cfg.frames_out, cfg.profile_out);

  // --- 3./4. CREA OSTACOLI E MASCHERA ------------------------------------
  // In 3D i due passi coincidono: le pareti sono le facce del dominio.
  const types::boundary_mask_t boundary_mask = build_boundary_mask(grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  //
  // Lo stesso writer va agganciato sia a `sim` sia a `solver`: il primo
  // notifica l'header con le dimensioni della griglia, il secondo i frame
  // delle norme di velocita', e sono due DataObservable distinti.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  LBMSimulation<DIM, D3Q27, COLLISION> simulation(
      grid_size, boundary_mask,
      CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, lid_velocity));
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D3Q27, COLLISION> solver(cfg.niters, cfg.nframes);
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

  return 0;
}
