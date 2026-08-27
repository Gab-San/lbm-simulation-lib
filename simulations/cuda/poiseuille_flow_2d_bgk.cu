// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"
#include "lbm/logging.hpp"

// QUILL LIB
#include "quill/LogMacros.h"

// C++ STD LIB
#include <memory>
#include <unordered_map>
#include <vector>

static constexpr unsigned short int DIM = 2;

template <unsigned short int dim> struct Config;
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

  /// Velocita' di riferimento, usata solo per calcolare nu/tau
  /// NON muove piu' nessuna parete in Poiseuille!!!
  const lbm::utils::Vector<double, 2> init_vel;

  /// Output path for frames
  const std::string out_frames;

  /// Output path for benchmark data
  const std::string out_data;

  const std::vector<lbm::CollisionDetection::CollisionArea<DIM>> obstacles;

  const std::unordered_map<unsigned int, uint8_t> obst_type_map;

  Config<2>(
      const lbm::types::DimPoint<2> grid_size_, const unsigned int c_iters,
      const unsigned int c_frames, const double c_reyn_num,
      const lbm::utils::Vector<double, 2> init_vel_,
      const std::string c_out_frames, const std::string c_out_data,
      const std::vector<lbm::CollisionDetection::CollisionArea<DIM>> obstacles_,
      const std::unordered_map<unsigned int, uint8_t> obst_type_map_)
      : grid_size(grid_size_), iters(c_iters), frames(c_frames),
        reyn_num(c_reyn_num), init_vel(init_vel_), out_frames(c_out_frames),
        out_data(c_out_data), obstacles(std::move(obstacles_)),
        obst_type_map(std::move(obst_type_map_)) {}
};

int main() {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  const Coordinate<DIM> ZERO(0, 0);
  const Coordinate<2> B129(0, 128);
  const Coordinate<2> C129(128, 128);
  const Coordinate<2> D129(128, 0);

  std::vector<Config<2>> configs{
      Config<2>(
          {129, 129}, /*iters*/ 100000, /*frames*/ 200, /*reyn*/ 100.0,
          /*init_vel*/ {0.1, 0}, "out/norms_poiseuille_cuda_129_100_01_bgk.bin",
          "out/data_poiseuille_cuda_129_100_01_bgk.bin",
          {
              CollisionDetection::CollisionArea(
                  ZERO,
                  {CollisionDetection::Segment(ZERO, D129),   // bottom (y=0)
                   CollisionDetection::Segment(B129, C129)}), // top (y=128)
              CollisionDetection::CollisionArea(              // LEFT WALL
                  ZERO,
                  {CollisionDetection::Segment(ZERO + Vector<int, DIM>(0, 1),
                                               B129 - Vector<int, DIM>(0, 1))}),
              CollisionDetection::CollisionArea( // RIGHT WALL
                  ZERO,
                  {CollisionDetection::Segment(C129 - Vector<int, DIM>(0, 1),
                                               D129 + Vector<int, DIM>(0, 1))}),
          },
          {{0, Solid::BB_RIGID_WALL}, // fixed top and bottom wall
           {1, Solid::PRESSURE_PERIODIC_INLET},
           {2, Solid::PRESSURE_PERIODIC_OUTLET}}), // right and left
                                                   // periodic bc
  };

  logging::setup_quill();
  quill::Logger *main_logger = logging::create_or_get_logger("main");

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  for (auto confidx = 0; confidx < configs.size(); confidx++) {
    const auto conf = configs[confidx];
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map] = conf;

    LOG_INFO(
        main_logger,
        "Simulation #{} Parameters:\n\tGrid dimensions: {}\n\tReynolds number: "
        "{}\n\tInitial Velocity: {}\n\tNumber of Iterations: {}\n\tNumber of "
        "frames: {}\n",
        confidx, grid_size, reyn, init_vel, iters, frames);

    types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<VtkWriter> writer =
        std::make_shared<VtkWriter>(conf.out_frames);

    CollisionParams<DIM, CollisionType> params(reyn, grid_size, init_vel);
    const double pout = 1;
    const double pin =
        pout + (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) *
                   8 * params.nu * params.init_vel.dx;
    Simulation simulation(grid_size, boundary_mask, params, pin, pout);

    simulation.attachListener(writer);

    CUDASolver<DIM, D2Q9, CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver);
    simulation.output(out_data.c_str(),
                      functional::extract_dx_profile_along_y_center);

    const double H = static_cast<double>(grid_size.y - 1);
    const auto exact_solution = analysis::PoiseuilleSolution2D(H, init_vel.dx);
    const double err_l2 =
        simulation.compute_error(analysis::NormType::L2, exact_solution);

    LOG_NOTICE(main_logger, "{} error: {}",
               analysis::to_string(analysis::NormType::L2), err_l2);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
