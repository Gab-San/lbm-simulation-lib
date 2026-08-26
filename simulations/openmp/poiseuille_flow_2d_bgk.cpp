// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/config/config-scanner.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"
#include "lbm/logging.hpp"

// QUILL LIB
#include "quill/LogMacros.h"

// C++ STD LIB
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
// Path assoluti ai dati di input, iniettati da CMake (vedi
// simulations/CMakeLists.txt): servono a non dipendere dalla working
// directory da cui si lancia l'eseguibile. I fallback valgono solo se si
// compila fuori da CMake.
#ifndef LBM_CONFIGS_DIR
#define LBM_CONFIGS_DIR "configs"
#endif

static constexpr unsigned short int DIM = 2;

/// Tipo di problema di cui si occupa questo main: i .toml con un
/// [problem].type diverso vengono ignorati.
static constexpr const char *PROBLEM_TYPE = "poiseuille";

/// L'operatore di collisione e' un parametro *template*, quindi va fissato a
/// compile time: questo binario esegue solo le configurazioni BGK, quelle
/// TRT le prende poiseuille_flow_2d_trt.
static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

/// Stessa logica per il backend.
static constexpr const char *BACKEND = "openmp";

/// Uso: poiseuille_flow_2d_bgk [cartella_configurazioni]
/// Senza argomenti usa la cartella configs/ dei sorgenti.
int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  const std::filesystem::path configs_dir =
      (argc > 1) ? argv[1] : LBM_CONFIGS_DIR;

  LOG_INFO(main_logger, "Reading configurations from '{}'",
           configs_dir.string());

  using Simulation = LBMSimulation<DIM, D2Q9, COLLISION>;

  unsigned int run_count = 0;

  config::run_matching_configs(
      configs_dir, PROBLEM_TYPE,
      [&](const config::SimulationConfig &cfg, const toml::table &) {
        if (cfg.collision != COLLISION || cfg.backend != BACKEND)
          return;

        const DimPoint<DIM> grid_size(cfg.nx, cfg.ny);
        const utils::Vector<double, DIM> init_vel(cfg.init_vel_x, 0.0);

        LOG_INFO(main_logger,
                 "Simulation '{}':\n\tGrid dimensions: {}\n\tReynolds number: "
                 "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: "
                 "{}\n\tNumber of frames: {}",
                 cfg.name, grid_size, cfg.reynolds, init_vel, cfg.niters,
                 cfg.nframes);

        // --- Geometria: resta qui, non nel TOML -------------------------
        // Poiseuille: pareti rigide sopra e sotto, ingresso e uscita a
        // pressione imposta sui lati (esclusi gli angoli, che appartengono
        // gia' alle pareti orizzontali). Tutto ricavato da nx/ny.
        const int x_max = static_cast<int>(cfg.nx) - 1;
        const int y_max = static_cast<int>(cfg.ny) - 1;

        const Coordinate<2> A(0, 0);
        const Coordinate<2> B(0, y_max);
        const Coordinate<2> C(x_max, y_max);
        const Coordinate<2> D(x_max, 0);

        const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles{
            CollisionDetection::CollisionArea(
                A, {CollisionDetection::Segment(A, D),   // bottom (y=0)
                    CollisionDetection::Segment(B, C)}), // top (y=y_max)
            CollisionDetection::CollisionArea(           // LEFT WALL
                A, {CollisionDetection::Segment(A + Vector<int, DIM>(0, 1),
                                                B - Vector<int, DIM>(0, 1))}),
            CollisionDetection::CollisionArea( // RIGHT WALL
                A, {CollisionDetection::Segment(C - Vector<int, DIM>(0, 1),
                                                D + Vector<int, DIM>(0, 1))}),
        };

        const std::unordered_map<unsigned int, uint8_t> obst_type_map{
            {0, Solid::BB_RIGID_WALL}, // fixed top and bottom wall
            {1, Solid::PRESSURE_PERIODIC_INLET},
            {2, Solid::PRESSURE_PERIODIC_OUTLET}}; // right and left
                                                   // periodic bc

        types::boundary_mask_t boundary_mask =
            Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles,
                                              grid_size);

        // --- Esecuzione --------------------------------------------------
        std::shared_ptr<AsyncBinaryWriter> writer =
            std::make_shared<AsyncBinaryWriter>(cfg.frames_file());

        const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size,
                                                     init_vel);
        // Salto di pressione che sostiene il flusso: ricavato dalla
        // soluzione di Poiseuille per il canale, non e' un parametro libero.
        const double pout = 1;
        const double pin =
            pout +
            (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) * 8 *
                params.nu * params.init_vel.dx;

        Simulation simulation(grid_size, boundary_mask, params, pin, pout);
        simulation.attachListener(writer);

        OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
        solver.attachListener(writer);

        simulation.solve(solver);

        simulation.output(cfg.profile_file().c_str(),
                          functional::extract_dx_profile_along_y_center);

        const double H = static_cast<double>(grid_size.y - 1);
        const auto exact_solution =
            analysis::PoiseuilleSolution2D(H, init_vel.dx);
        const double err_l2 =
            simulation.compute_error(analysis::NormType::L2, exact_solution);

        LOG_NOTICE(main_logger, "{} error: {}",
                   analysis::to_string(analysis::NormType::L2), err_l2);

        simulation.detachListener(writer);
        solver.detachListener(writer);
        ++run_count;
      });

  if (run_count == 0) {
    LOG_WARNING(main_logger,
                "Nessuna configurazione eseguita: in '{}' non c'e' nessun "
                "file .toml con [problem].type = \"{}\", [physics].collision = "
                "\"{}\" e [solver].backend = \"{}\"",
                configs_dir.string(), PROBLEM_TYPE,
                collision_model_to_string(COLLISION), BACKEND);
    return 1;
  }

  LOG_INFO(main_logger, "Number of Simulations: {}", run_count);

  return 0;
}
