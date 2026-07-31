// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/problems/problem_2d.hpp"
#include "lbm-sim/solver/solver-2d.hpp"
#include "lbm-sim/structure.hpp"

// COLLISION DETECTION LIB
#include "collision-detection/collision-area.hpp"
#include "collision-detection/core/types.hpp"

// C++ STD LIB
#include <memory>

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

  const std::vector<CollisionDetection::CollisionArea<2>> obstacles;

  Config<2>(
      const std::size_t c_x, const std::size_t c_y, const unsigned int c_iters,
      const unsigned int c_frames, const double c_reyn_num,
      const double c_init_vel, const std::string c_out_frames,
      const std::string c_out_data,
      const std::vector<CollisionDetection::CollisionArea<2>> &&obstacles_)
      : grid_num_cells_x(c_x), grid_num_cells_y(c_y), iters(c_iters),
        frames(c_frames), reyn_num(c_reyn_num), init_vel(c_init_vel),
        out_frames(c_out_frames), out_data(c_out_data), obstacles(obstacles_) {}
};

int main() {
  using namespace lbm;
  using CollisionDetection::types::Coordinate;
  using CollisionDetection::types::DimPoint;

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(128, 128);
  const Coordinate<2> D(128, 0);

  std::vector<Config<DIM>> configs{
      Config<DIM>(129, 129, /*iters*/ 10000, /*frames*/ 100,
                  /*reyn*/ 100.0, /*init_vel*/ 0.1,
                  "out/norms_129_100_01_bgk.bin", "out/data_129_100_01_bgk.bin",
                  {CollisionDetection::CollisionArea(
                       A, {CollisionDetection::Segment(A, B),
                           CollisionDetection::Segment(A, D),
                           CollisionDetection::Segment(D, C)}),
                   CollisionDetection::CollisionArea(
                       A, {CollisionDetection::Segment(B, C)})})};

  constexpr auto CollisionType = CollisionModel::BGK;

  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[size_x, size_y, iters, frames, reyn, init_vel, out_frames,
                 out_data, obstacles] = conf;

    const Structure<DIM> strt(obstacles, 1);
    strt.obstacles[0].getPerimeter();
    strt.obstacles[1].getPerimeter();

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(out_frames);

    Simulation simulation(
        DimPoint<DIM>(size_x, size_y),
        Params<DIM, CollisionType>(reyn, {size_x, size_y}, {init_vel, 0.0}));

    simulation.attachListener(writer);

    MPISolver2D<CollisionType> solver(iters, frames, strt);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/, problem);

    simulation.output(out_data.c_str());

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
