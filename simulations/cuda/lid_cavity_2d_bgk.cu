#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"
#include "lbm/logging.hpp"

// QUILL LIB
#include "quill/LogMacros.h"

// C++ STD LIB

static constexpr unsigned short int DIM = 2;

template <unsigned short int dim> struct Config;

/**
 * \brief This struct represents a configuration.
 * A configuration is an instance of a simulation.
 *
 * The parameters of a lid cavity simulation are:
 * - grid_size (num_cells_x, num_cells_y)
 * - reyn_num: reynold number
 * - u_lid: initial velocity of the lid
 * - coll_op: collision operator
 * - num_steps: number of iteration steps
 * - num_frames: number of frames to save
 *
 *
 * Frames contain the information about
 * the norm of the velocity at a step t.
 *
 *
 * Jobs are divided from one another by
 * whitelines.
 * There is no special string of characters to divide jobs,
 * but all of the parameters of a job must be configured.
 *
 * \note TRT operator is not supported
 *
 */
template <> struct Config<2> {

  const lbm::types::DimPoint<2> grid_size;

  /// Number of iteration steps
  const unsigned int iters;

  /// Number of frames
  ///
  /// Frames contain the information about
  /// the norm of the velocity at a step t.
  const unsigned int frames;

  /// Reynold number
  const double reyn_num;

  /// Initial velocity of the fluid
  const lbm::utils::Vector<double, 2> init_vel;

  /// Output path for frames
  const std::string out_frames;

  /// Output path for benchmark data
  const std::string out_data;

  const std::vector<lbm::CollisionDetection::CollisionArea<DIM>> obstacles;

  /// Tabella laterale: id ostacolo -> {tipo di BC, velocita' di parete}.
  const std::vector<lbm::Solid::ObstacleData<DIM>> obstacle_data;

  /// BC delle facce del dominio: le pareti non sono piu' ostacoli.
  const lbm::Solid::DomainBC<DIM> domain_bc;

  Config<2>(
      const lbm::types::DimPoint<2> grid_size_, const unsigned int c_iters,
      const unsigned int c_frames, const double c_reyn_num,
      const lbm::utils::Vector<double, 2> init_vel_,
      const std::string c_out_frames, const std::string c_out_data,
      const std::vector<lbm::CollisionDetection::CollisionArea<DIM>> obstacles_,
      const std::vector<lbm::Solid::ObstacleData<DIM>> obstacle_data_,
      const lbm::Solid::DomainBC<DIM> domain_bc_)
      : grid_size(grid_size_), iters(c_iters), frames(c_frames),
        reyn_num(c_reyn_num), init_vel(init_vel_), out_frames(c_out_frames),
        out_data(c_out_data), obstacles(std::move(obstacles_)),
        obstacle_data(std::move(obstacle_data_)), domain_bc(domain_bc_) {}
};

/// Lid-driven cavity: tre pareti rigide + il lid mobile in alto.
static lbm::Solid::DomainBC<DIM> make_cavity_bc() {
  lbm::Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = lbm::Solid::BB_RIGID_WALL;   // x = 0
  dbc.high(0) = lbm::Solid::BB_RIGID_WALL;  // x = nx-1
  dbc.low(1) = lbm::Solid::BB_RIGID_WALL;   // y = 0
  dbc.high(1) = lbm::Solid::BB_MOVING_WALL; // il lid, y = ny-1
  return dbc;
}

int main() {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;

  const Coordinate<2> ZERO(0, 0);

  const Coordinate<2> B200(0, 199);
  const Coordinate<2> C200(199, 199);
  const Coordinate<2> D200(199, 0);

  std::vector<Config<DIM>> configs{
      Config<DIM>({129, 129}, /*iters*/ 10000, /*frames*/ 100,
                  /*reyn*/ 100.0, /*init_vel*/ {0.1, 0},
                  "out/norms_lid_cavity_cuda_129_100_01_bgk.bin",
                  "out/data_lid_cavity_cuda_129_100_01_bgk.bin", {}, {},
                  make_cavity_bc()),

      Config<DIM>({200, 200}, /*iters*/ 30000, /*frames*/ 100,
                  /*reyn*/ 1000.0, /*init_vel*/ {0.1, 0},
                  "out/norms_lid_cavity_cuda_200_1000_01_bgk.bin",
                  "out/data_lid_cavity_cuda_200_1000_01_bgk.bin", {}, {},
                  make_cavity_bc()),

#if 0
      Config<DIM>(
          {2000, 2000}, /*iters*/ 50000, /*frames*/ 200, /*reyn*/ 7500,
          /*init_vel*/ {0.2, 0},
          "out/norms_lid_cavity_cuda_2000_7500_02_bgk.bin",
          "out/data_lid_cavity_cuda_2000_7500_02_bgk.bin",
          {}, {}, make_cavity_bc()),

      Config<DIM>(
          {5000, 5000}, /*iters*/ 100000, /*frames*/ 250,
          /*reyn*/ 50000.0, /*init_vel*/ {0.1, 0},
          "out/norms_lid_cavity_cuda_5000_50000_01_bgk.bin",
          "out/data_lid_cavity_cuda_5000_50000_01_bgk.bin",
          {}, {}, make_cavity_bc()),
#endif
  };

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  std::string path_to_benchmark("benchmarks/ghia/");

  LOG_INFO(main_logger, "Number of Simulations: {}", configs.size());

  for (std::size_t confidx = 0; confidx < configs.size(); confidx++) {
    const auto conf = configs[confidx];
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obstacle_data, domain_bc] = conf;

    LOG_INFO(
        main_logger,
        "Simulation #{} Parameters:\n\tGrid dimensions: {}\n\tReynolds number: "
        "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber of "
        "frames: {}\n",
        confidx, grid_size, reyn, init_vel, iters, frames);

    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>(obstacles, grid_size);

    std::shared_ptr<VtkWriter> writer = std::make_shared<VtkWriter>(out_frames);

    Simulation simulation(
        grid_size, std::move(solid_mask), obstacle_data, domain_bc,
        CollisionParams<DIM, CollisionType>(reyn, grid_size, init_vel));

    simulation.attachListener(writer);

    CUDASolver<DIM, D2Q9, CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/);

    simulation.output(out_data.c_str(),
                      functional::extract_dy_profile_along_x_center);

    // Confronto con Ghia et al. (1982): Norma scelta qui: L2.
    const auto ghia_y = simulation.compute_ghia_error(
        path_to_benchmark + "data_y_" + formatting::format_reyn(reyn) + ".txt");

    LOG_NOTICE(main_logger, "Ghia ({}) | uy(x/2): rel={} abs={}",
               analysis::to_string(analysis::NormType::L2), ghia_y.relative,
               ghia_y.absolute);

    const auto ghia_x = simulation.compute_ghia_error(
        path_to_benchmark + "data_x_" + formatting::format_reyn(reyn) + ".txt");

    LOG_NOTICE(main_logger, "Ghia ({}) | ux(y/2): rel={} abs={}",
               analysis::to_string(analysis::NormType::L2), ghia_x.relative,
               ghia_x.absolute);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
