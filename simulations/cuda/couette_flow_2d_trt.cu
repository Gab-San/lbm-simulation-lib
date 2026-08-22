// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/cuda-solver.cuh"

#include "lbm-sim/collision-detection/collision-area.hpp"

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
  using utils::Vector;

  const Coordinate<2> ZERO(0, 0);
  const Coordinate<2> B129(0, 128);
  const Coordinate<2> C129(128, 128);
  const Coordinate<2> D129(128, 0);

  std::vector<Config<2>> configs{
      Config<2>(
          {129, 129}, /*iters*/ 130000, /*frames*/ 300, /*reyn*/ 100.0,
          /*init_vel*/ {0.1, 0}, "out/norms_couette_openmp_129_100_01.bin",
          "out/data_couette_cuda_129_100_01.bin",
          {
              CollisionDetection::CollisionArea(
                  ZERO,
                  {CollisionDetection::Segment(ZERO, D129)}), // bottom (y=0)
              CollisionDetection::CollisionArea(
                  ZERO,
                  {CollisionDetection::Segment(B129, C129)}), // top (y=128)
              CollisionDetection::CollisionArea(
                  ZERO,
                  {CollisionDetection::Segment(ZERO + Vector<int, DIM>(0, 1),
                                               B129 - Vector<int, DIM>(0, 1)),
                   CollisionDetection::Segment(C129 - Vector<int, DIM>(0, 1),
                                               D129 + Vector<int, DIM>(0, 1))}),
          },
          {{0, Solid::BB_RIGID_WALL},
           {1, Solid::BB_MOVING_WALL},
           {2, Solid::PERIODIC}}),
  };

  constexpr auto CollisionType = CollisionModel::TRT;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;
  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map] = conf;
    types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(conf.out_frames);

    Simulation simulation(
        grid_size, boundary_mask,
        CollisionParams<DIM, CollisionType>(reyn, grid_size, init_vel));

    simulation.attachListener(writer);

    CUDASolver2D<CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver, problem);
    simulation.output(out_data.c_str());

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
