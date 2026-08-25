#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/problems/problem_2d.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"

// C++ STD LIB
#include <unordered_map>

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

  const Coordinate<2> ZERO(0, 0);

  const Coordinate<2> B200(0, 199);
  const Coordinate<2> C200(199, 199);
  const Coordinate<2> D200(199, 0);

  std::vector<Config<DIM>> configs{
      Config<DIM>(
          {129, 129}, /*iters*/ 10000, /*frames*/ 100,
          /*reyn*/ 100.0, /*init_vel*/ {0.1, 0},
          "out/norms_lid_cavity_cuda_129_100_01_bgk.bin",
          "out/data_lid_cavity_cuda_129_100_bgk.bin",
          {CollisionDetection::CollisionArea(
               ZERO,
               {CollisionDetection::Segment(ZERO, Coordinate<DIM>(0, 128)),
                CollisionDetection::Segment(ZERO, Coordinate<DIM>(128, 0)),
                CollisionDetection::Segment(Coordinate<DIM>(128, 0),
                                            Coordinate<DIM>(128, 128))}),
           CollisionDetection::CollisionArea(
               ZERO, {CollisionDetection::Segment(Coordinate<DIM>(0, 128),
                                                  Coordinate<DIM>(128, 128))})},
          {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),

      Config<DIM>({200, 200}, /*iters*/ 30000, /*frames*/ 100,
                  /*reyn*/ 1000.0, /*init_vel*/ {0.1, 0},
                  "out/norms_lid_cavity_cuda_200_1000_01_bgk.bin",
                  "out/data_lid_cavity_cuda_200_1000_01_bgk.bin",
                  {CollisionDetection::CollisionArea(
                       ZERO, {CollisionDetection::Segment(ZERO, B200),
                              CollisionDetection::Segment(ZERO, D200),
                              CollisionDetection::Segment(D200, C200)}),
                   CollisionDetection::CollisionArea(
                       ZERO, {CollisionDetection::Segment(B200, C200)})},
                  {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),

#if 0
      Config<DIM>(
          {2000, 2000}, /*iters*/ 50000, /*frames*/ 200, /*reyn*/ 7500,
          /*init_vel*/ {0.2, 0},
          "out/norms_lid_cavity_cuda_2000_7500_02_bgk.bin",
          "out/data_lid_cavity_cuda_2000_7500_02_bgk.bin",
          {CollisionDetection::CollisionArea(
               ZERO,
               {CollisionDetection::Segment(ZERO, Coordinate<DIM>(0, 1999)),
                CollisionDetection::Segment(ZERO, Coordinate<DIM>(1999, 0)),
                CollisionDetection::Segment(Coordinate<DIM>(1999, 0),
                                            Coordinate<DIM>(1999, 1999))}),
           CollisionDetection::CollisionArea(
               ZERO,
               {CollisionDetection::Segment(Coordinate<DIM>(0, 1999),
                                            Coordinate<DIM>(1999, 1999))})},
          {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),

      Config<DIM>(
          {5000, 5000}, /*iters*/ 100000, /*frames*/ 250,
          /*reyn*/ 50000.0, /*init_vel*/ {0.1, 0},
          "out/norms_lid_cavity_cuda_5000_50000_01_bgk.bin",
          "out/data_lid_cavity_cuda_5000_50000_01_bgk.bin",
          {CollisionDetection::CollisionArea(
               ZERO,
               {CollisionDetection::Segment(ZERO, Coordinate<DIM>(0, 4999)),
                CollisionDetection::Segment(ZERO, Coordinate<DIM>(4999, 0)),
                CollisionDetection::Segment(Coordinate<DIM>(4999, 0),
                                            Coordinate<DIM>(4999, 4999))}),
           CollisionDetection::CollisionArea(
               ZERO,
               {CollisionDetection::Segment(Coordinate<DIM>(0, 4999),
                                            Coordinate<DIM>(4999, 4999))})},
          {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),
#endif
  };

  constexpr auto CollisionType = CollisionModel::BGK;

  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map] = conf;

    types::boundary_mask_t obstacle_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(out_frames);

    Simulation simulation(
        grid_size, obstacle_mask,
        CollisionParams<DIM, CollisionType>(reyn, grid_size, init_vel));

    simulation.attachListener(writer);

    CUDASolver2D<CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/, problem);

    simulation.output(out_data.c_str(),
                      functional::extract_dy_profile_along_x_center);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
