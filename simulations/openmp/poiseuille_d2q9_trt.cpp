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
#include "lbm/format/csv-writer.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"

// C++ STD LIB
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr lbm::types::dim_t DIM = 2;
constexpr auto COLLISION = lbm::CollisionModel::TRT;
constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;

  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

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

  LBM_LOG_INFO(
      main_logger,
      "Simulation:\n\tGrid dimensions: {}\n\tReynolds number: "
      "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
      "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
      grid_size, cfg.reynolds, cfg.u0, cfg.niters, cfg.nframes, cfg.frames_out,
      cfg.profile_out);

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
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({}, grid_size);

  // --- 5. LANCIA SIMULAZIONE ---------------------------------------------
  // frames_out e' la CARTELLA; il basename dei file lo da' il nome della
  // configurazione, cosi' run diversi nella stessa cartella non si
  // sovrascrivono a vicenda.
  std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

  const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size, cfg.u0);
  // Salto di pressione che sostiene il flusso: ricavato dalla soluzione di
  // Poiseuille per il canale, non e' un parametro libero.
  const double pout = 1;
  const double pin =
      pout +
      numbers::invcs_2 *
          (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) * 8 *
          params.nu * params.init_vel.dx;

  LBMSimulation<DIM, D2Q9, COLLISION> simulation(
      grid_size, std::move(solid_mask), {}, dbc, params, pin, pout);
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
                    functional::extract_dx_profile_along_y_center);

  // --- 7. CALCOLO DELL'ERRORE --------------------------------------------
  // H = altezza canale (parete inferiore a y=0, superiore a y=ny-1);
  // Umax = velocita' di riferimento (parete mobile per Couette).
  const double H = static_cast<double>(grid_size.y - 1);
  const auto exact_solution = analysis::PoiseuilleSolution2D(H, cfg.u0.dx);
  const auto error_metrics = simulation.compute_error_analysis(
      analysis::NormType::L2, exact_solution, cfg.u0.dx);
  const double err_l2 = error_metrics.absolute;

  LBM_LOG_NOTICE(main_logger, "{} error: {}",
                 analysis::to_string(analysis::NormType::L2), err_l2);

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
      "Poiseuille validation summary\n"
      "  Grid:                 {} x {}\n"
      "  Reynolds:             {}\n"
      "  Collision:            {}\n"
      "  Iterations:           {}\n\n"
      "  Error vs analytical profile:\n"
      "    absolute L2:        {:.8e}\n"
      "    relative L2:        {:.4f} %\n"
      "    RMSE:               {:.8e}\n"
      "    RMSE / U_max:       {:.4f} %\n"
      "    Linf:               {:.8e}\n"
      "    Linf / U_max:       {:.4f} %\n\n"
      "  Runtime:              {:.3f} s\n"
      "  Threads:              {}\n"
      "  Performance:          {:.3f} MLUPS",
      grid_size.x, grid_size.y, cfg.reynolds,
      collision_model_to_string(COLLISION), cfg.niters, error_metrics.absolute,
      100.0 * error_metrics.relative, error_metrics.rmse,
      100.0 * error_metrics.rmse_normalized, error_metrics.linf,
      100.0 * error_metrics.linf_normalized, runtime_s, n_threads, mlups);

  format::CsvWriter<analysis::DetailedErrorAnalysisSchema> error_writer(
      "out/error_" + cfg.name + "_" + format::format_reyn(cfg.reynolds) +
          ".csv",
      true);

  error_writer.append_row(
      "ux", grid_size.x, grid_size.y, cfg.reynolds,
      collision_model_to_string(COLLISION), cfg.niters,
      analysis::to_string(analysis::NormType::L2), error_metrics.relative,
      error_metrics.absolute, 100.0 * error_metrics.relative,
      error_metrics.rmse, 100.0 * error_metrics.rmse_normalized,
      error_metrics.linf, 100.0 * error_metrics.linf_normalized, runtime_s,
      n_threads, mlups);

  error_writer.flush();
  error_writer.close();

  simulation.detachListener(writer);
  solver.detachListener(writer);
           #ifdef LBM_PROFILING
  lbm::profiling::dump_csv(cfg.profile_out);  
  lbm::profiling::reset();                    
#endif
}
  return 0;
}
