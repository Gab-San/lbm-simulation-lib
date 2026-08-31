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

    const utils::Vector<double, DIM> init_vel = cfg.u0;

    LBM_LOG_INFO(
        main_logger,
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
    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

    LBMSimulation<DIM, D2Q9, COLLISION> simulation(
        grid_size, std::move(solid_mask), {}, dbc,
        CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, init_vel));
    simulation.attachListener(writer);

    OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    const auto solve_start = std::chrono::steady_clock::now();
    simulation.solve(solver);
    const auto solve_end = std::chrono::steady_clock::now();
    const double runtime_s =
        std::chrono::duration<double>(solve_end - solve_start).count();

    // --- 6. OUTPUT ---------------------------------------------------------
    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dy_profile_along_x_center);

    // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
    // Confronto con Ghia et al. (1982). Norma scelta qui: L2.
    /* const std::string path_to_benchmark =
        std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

    const auto ghia_y = simulation.compute_ghia_error(
        path_to_benchmark + "data_y_" + format::format_reyn(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_y.relative,
                   ghia_y.absolute);

    const auto ghia_x = simulation.compute_ghia_error(
        path_to_benchmark + "data_x_" + format::format_reyn(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_x.relative,
                   ghia_x.absolute);

    simulation.detachListener(writer);
    solver.detachListener(writer);

#ifdef LBM_PROFILING
    lbm::profiling::dump_csv(cfg.profile_out);
    lbm::profiling::reset();

#endif
    */
   const std::string path_to_benchmark =
        std::string(LBM_BENCHMARKS_DIR) + "/ghia/";

    format::CsvWriter<analysis::DetailedErrorAnalysisSchema> error_writer(
        "out/error_" + cfg.name + "_" + format::format_reyn(cfg.reynolds) +
            ".csv",
        true);

    const auto ghia_y = simulation.compute_ghia_error(
        path_to_benchmark + "data_y_" + format::format_reyn(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_y.relative,
                   ghia_y.absolute);

    const auto ghia_x = simulation.compute_ghia_error(
        path_to_benchmark + "data_x_" + format::format_reyn(cfg.reynolds) +
        ".txt");

    LBM_LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
                   analysis::to_string(analysis::NormType::L2), ghia_x.relative,
                   ghia_x.absolute);

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads() >= omp_get_num_procs()
                    ? omp_get_num_procs()
                    : omp_get_max_threads();
#endif

    const double mlups =
        runtime_s > 0.0
            ? (static_cast<double>(grid_size.x) *
               static_cast<double>(grid_size.y) *
               static_cast<double>(cfg.niters)) /
                  (runtime_s * 1.0e6)
            : 0.0;

    LBM_LOG_NOTICE(
        main_logger,
        "Lid-driven cavity validation summary\n"
        "  Grid:                 {} x {}\n"
        "  Reynolds:             {}\n"
        "  Collision:            {}\n"
        "  Iterations:           {}\n\n"
        "  Ghia ux(y/2):\n"
        "    relative L2:        {:.4f} %\n"
        "    RMSE / U_lid:       {:.4f} %\n"
        "    Linf / U_lid:       {:.4f} %\n\n"
        "  Ghia uy(x/2):\n"
        "    relative L2:        {:.4f} %\n"
        "    RMSE / U_lid:       {:.4f} %\n"
        "    Linf / U_lid:       {:.4f} %\n\n"
        "  Runtime:              {:.3f} s\n"
        "  Threads:              {}\n"
        "  Performance:          {:.3f} MLUPS",
        grid_size.x, grid_size.y, cfg.reynolds,
        collision_model_to_string(COLLISION), cfg.niters,
        100.0 * ghia_x.relative, 100.0 * ghia_x.rmse_normalized,
        100.0 * ghia_x.linf_normalized, 100.0 * ghia_y.relative,
        100.0 * ghia_y.rmse_normalized, 100.0 * ghia_y.linf_normalized,
        runtime_s, n_threads, mlups);

    error_writer.append_row("uy", grid_size.x, grid_size.y, cfg.reynolds,
                            collision_model_to_string(COLLISION), cfg.niters,
                            analysis::to_string(analysis::NormType::L2),
                            ghia_y.relative, ghia_y.absolute,
                            100.0 * ghia_y.relative, ghia_y.rmse,
                            100.0 * ghia_y.rmse_normalized, ghia_y.linf,
                            100.0 * ghia_y.linf_normalized, runtime_s,
                            n_threads, mlups);

    error_writer.append_row("ux", grid_size.x, grid_size.y, cfg.reynolds,
                            collision_model_to_string(COLLISION), cfg.niters,
                            analysis::to_string(analysis::NormType::L2),
                            ghia_x.relative, ghia_x.absolute,
                            100.0 * ghia_x.relative, ghia_x.rmse,
                            100.0 * ghia_x.rmse_normalized, ghia_x.linf,
                            100.0 * ghia_x.linf_normalized, runtime_s,
                            n_threads, mlups);

    error_writer.flush();
    error_writer.close();
  }
  return 0;
}
