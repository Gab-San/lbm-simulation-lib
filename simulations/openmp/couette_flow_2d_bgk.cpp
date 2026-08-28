// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/config/config-parser.hpp"
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

static constexpr lbm::types::dim_t DIM = 2;

static constexpr lbm::CollisionModel COLLISION = lbm::CollisionModel::BGK;
static constexpr lbm::ExecutionBackend BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  // --- 1. LEGGI CONFIGURAZIONI --------------------------------------------
  if (argc < 2) {
    config::print_usage(argv[0]);
    return 1;
  }

  std::vector<config::SimulationConfig<DIM>> configs;
  try {
    configs = config::parse_config<DIM>(argv[1]);
    for (auto &cfg : configs)
      config::ensure_compatible(cfg, COLLISION, BACKEND);
  } catch (const config::ConfigError &err) {
    std::cerr << "Errore di configurazione: " << err.what() << "\n";
    return 1;
  }

  // --- 2. ISTANZIA LOGGER --------------------------------------------------
  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  // --- 3. ESEGUI UNA SIMULAZIONE PER OGNI CONFIG ---------------------------
  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);

    LOG_INFO(main_logger,
             "Simulation:\n\tGrid dimensions: {}\n\tReynolds number: "
             "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber "
             "of frames: {}\n\tFrames output: {}\n\tProfile output: {}",
             grid_size, cfg.reynolds, cfg.u0, cfg.niters, cfg.nframes,
             cfg.frames_out, cfg.profile_out);

    const int x_max = static_cast<int>(grid_size.x) - 1;
    const int y_max = static_cast<int>(grid_size.y) - 1;

    const Coordinate<2> A(0, 0);
    const Coordinate<2> B(0, y_max);
    const Coordinate<2> C(x_max, y_max);
    const Coordinate<2> D(x_max, 0);

    const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles{
        CollisionDetection::CollisionArea(A,
                                          {CollisionDetection::Segment(A, D)}),
        CollisionDetection::CollisionArea(A,
                                          {CollisionDetection::Segment(B, C)}),
        CollisionDetection::CollisionArea(
            A, {CollisionDetection::Segment(A + Vector<int, DIM>(0, 1),
                                            B - Vector<int, DIM>(0, 1)),
                CollisionDetection::Segment(C - Vector<int, DIM>(0, 1),
                                            D + Vector<int, DIM>(0, 1))}),
    };

    const std::unordered_map<unsigned int, uint8_t> obst_type_map{
        {0, Solid::BB_RIGID_WALL},
        {1, Solid::BB_MOVING_WALL},
        {2, Solid::PERIODIC}};

    types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<VtkWriter> writer =
        std::make_shared<VtkWriter>(cfg.frames_out);

    LBMSimulation<DIM, D2Q9, COLLISION> simulation(
        grid_size, boundary_mask,
        CollisionParams<DIM, COLLISION>(cfg.reynolds, grid_size, cfg.u0));
    simulation.attachListener(writer);

    OpenMPSolver<DIM, D2Q9, COLLISION> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver);

    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dx_profile_along_y_center);

    const double H = static_cast<double>(grid_size.y - 1);
    const auto exact_solution = analysis::CouetteSolution2D(H, cfg.u0.dx);
    const double err_l2 =
        simulation.compute_error(analysis::NormType::L2, exact_solution);

    LOG_NOTICE(main_logger, "{} error: {}",
               analysis::to_string(analysis::NormType::L2), err_l2);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
