// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/problems/problem_2d.hpp"
#include "lbm-sim/solver/solver-2d.hpp"

// COLLISION DETECTION LIB
#include "collision-detection/collision-area.hpp"
#include "collision-detection/core/types.hpp"

// C++ STD LIB
#include <memory>
#include <unordered_map>

static constexpr unsigned int DIM = 2;

template <unsigned int dim> struct Config;

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
  /// Number of cells along the x axis
  const std::size_t grid_num_cells_x;

  /// Number of cells along the y axis
  const std::size_t grid_num_cells_y;

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
  const double init_vel;

  /// Output path for frames
  const std::string out_frames;

  /// Output path for benchmark data
  const std::string out_data;

  const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles;

  const std::unordered_map<unsigned int, uint8_t> strt;

  Config<2>(
      const std::size_t c_x, const std::size_t c_y, const unsigned int c_iters,
      const unsigned int c_frames, const double c_reyn_num,
      const double c_init_vel, const std::string c_out_frames,
      const std::string c_out_data,
      const std::vector<CollisionDetection::CollisionArea<DIM>> &obstacles_,
      const std::unordered_map<unsigned int, uint8_t> &strt_)
      : grid_num_cells_x(c_x), grid_num_cells_y(c_y), iters(c_iters),
        frames(c_frames), reyn_num(c_reyn_num), init_vel(c_init_vel),
        out_frames(c_out_frames), out_data(c_out_data), obstacles(obstacles_),
        strt(strt_) {}
};

int main() {
  using namespace lbm;
  using CollisionDetection::types::Coordinate;
  using CollisionDetection::types::DimPoint;

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(128, 128);
  const Coordinate<2> D(128, 0);

  const Coordinate<2> A2(0, 0);
  const Coordinate<2> B2(0, 199);
  const Coordinate<2> C2(199, 199);
  const Coordinate<2> D2(199, 0);

  const Coordinate<2> B_ = B2 - Coordinate<2>(0, 1);
  const Coordinate<2> C_ = C2 - Coordinate<2>(0, 1);

  std::vector<Config<DIM>> configs{
      Config<DIM>(129, 129, /*iters*/ 10000, /*frames*/ 100,
                  /*reyn*/ 100.0, /*init_vel*/ 0.1,
                  "out/norms_129_100_01_bgk.bin", "out/data_129_100_01_bgk.bin",
                  {CollisionDetection::CollisionArea(
                       A, {CollisionDetection::Segment(A, B),
                           CollisionDetection::Segment(A, D),
                           CollisionDetection::Segment(D, C)}),
                   CollisionDetection::CollisionArea(
                       A, {CollisionDetection::Segment(B, C)})},
                  {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),

      Config<DIM>(200, 200, /*iters*/ 30000, /*frames*/ 100,
                  /*reyn*/ 1000.0, /*init_vel*/ 0.1,
                  "out/norms_200_1000_01_bgk.bin",
                  "out/data_200_1000_01_bgk.bin",
                  {CollisionDetection::CollisionArea(
                       A2, {CollisionDetection::Segment(A2, B2),
                            CollisionDetection::Segment(A2, D2),
                            CollisionDetection::Segment(D2, C2)}),
                   CollisionDetection::CollisionArea(
                       A2, {CollisionDetection::Segment(B2, C2)})},
                  {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),
  };

  constexpr auto CollisionType = CollisionModel::BGK;

  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[size_x, size_y, iters, frames, reyn, init_vel, out_frames,
                 out_data, obstacles, strt] = conf;

    Solid::types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(strt, obstacles, {size_x, size_y});

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(out_frames);

    Simulation simulation(
        DimPoint<DIM>(size_x, size_y),
        Params<DIM, CollisionType>(reyn, {size_x, size_y}, {init_vel, 0.0}));

    simulation.attachListener(writer);

    MPISolver2D<CollisionType> solver(iters, frames, boundary_mask);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/, problem);

    simulation.output(out_data.c_str());

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
