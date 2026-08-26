#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/config/config-scanner.hpp"
#include "lbm-sim/core/vector.hpp"
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
#ifndef LBM_BENCHMARKS_DIR
#define LBM_BENCHMARKS_DIR "benchmarks"
#endif

static constexpr lbm::types::dim_t DIM = 2;

/// Tipo di problema di cui si occupa questo main: i .toml con un
/// [problem].type diverso vengono ignorati.
static constexpr const char *PROBLEM_TYPE = "lid_cavity";

/// L'operatore di collisione e' un parametro *template*, quindi va fissato a
/// compile time: questo binario esegue solo le configurazioni BGK, quelle
/// TRT le prende lid_cavity_2d_trt.
static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;

/// Stessa logica per il backend: le configurazioni CUDA hanno griglie e
/// numeri di iterazioni pensati per la GPU, non vanno eseguite qui.
static constexpr const char *BACKEND = "openmp";

/// Uso: lid_cavity_2d_bgk [cartella_configurazioni]
/// Senza argomenti usa la cartella configs/ dei sorgenti.
int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;

  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  const std::filesystem::path configs_dir =
      (argc > 1) ? argv[1] : LBM_CONFIGS_DIR;
  const std::string path_to_benchmark =
      std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

  LOG_INFO(main_logger, "Reading configurations from '{}'",
           configs_dir.string());

  using Simulation = LBMSimulation<DIM, D2Q9, COLLISION>;

  unsigned int run_count = 0;

  config::run_matching_configs(
      configs_dir, PROBLEM_TYPE,
      [&](const config::SimulationConfig &cfg, const toml::table &) {
        // Un .toml di tipo lid_cavity ma con un altro operatore di
        // collisione, o destinato a un altro backend, e' di un altro binario.
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
        // Per la lid cavity gli "ostacoli" sono i quattro lati del dominio,
        // quindi si ricavano da nx/ny: tre pareti rigide + il lid mobile
        // in alto. E' questa la parte che dipende dal tipo di problema.
        const int x_max = static_cast<int>(cfg.nx) - 1;
        const int y_max = static_cast<int>(cfg.ny) - 1;

        const Coordinate<2> A(0, 0);
        const Coordinate<2> B(0, y_max);
        const Coordinate<2> C(x_max, y_max);
        const Coordinate<2> D(x_max, 0);

        const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles{
            CollisionDetection::CollisionArea(
                A, {CollisionDetection::Segment(A, B),
                    CollisionDetection::Segment(A, D),
                    CollisionDetection::Segment(D, C)}),
            CollisionDetection::CollisionArea(
                A, {CollisionDetection::Segment(B, C)})};

        const std::unordered_map<unsigned int, uint8_t> obst_type_map{
            {0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}};

        types::boundary_mask_t obstacle_mask =
            Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles,
                                              grid_size);

        // --- Esecuzione --------------------------------------------------
        std::shared_ptr<AsyncBinaryWriter> writer =
            std::make_shared<AsyncBinaryWriter>(cfg.frames_file());

        Simulation simulation(
            grid_size, obstacle_mask,
            CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, init_vel));
        simulation.attachListener(writer);

        OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
        solver.attachListener(writer);

        simulation.solve(solver /*, preconditioner*/);

        simulation.output(cfg.profile_file().c_str(),
                          functional::extract_dy_profile_along_x_center);

        // Confronto con Ghia et al. (1982): Norma scelta qui: L2.
        const auto ghia_y = simulation.compute_ghia_error(
            path_to_benchmark + "data_y_" +
            formatting::format_reyn(cfg.reynolds) + ".txt");

        LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_y.relative,
                   ghia_y.absolute);

        const auto ghia_x = simulation.compute_ghia_error(
            path_to_benchmark + "data_x_" +
            formatting::format_reyn(cfg.reynolds) + ".txt");

        LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_x.relative,
                   ghia_x.absolute);

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
