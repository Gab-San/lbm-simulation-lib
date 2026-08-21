#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/solver-2d.hpp"

#include "lbm-sim/collision-detection/collision-area.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

// C++ STD LIB
#include <memory>
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

  Config(
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

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(128, 128);
  const Coordinate<2> D(128, 0);

  const Coordinate<2> A2(0, 0);
  const Coordinate<2> B2(0, 199);
  const Coordinate<2> C2(199, 199);
  const Coordinate<2> D2(199, 0);

  std::vector<Config<DIM>> configs{
      Config<DIM>({129, 129}, /*iters*/ 10000, /*frames*/ 100,
                  /*reyn*/ 100.0, /*init_vel*/ {0.1, 0},
                  "out/norms_lid_cavity_openmp_129_100_01_bgk.bin",
                  "out/data_lid_cavity_openmp_129_100_01_bgk.bin",
                  {CollisionDetection::CollisionArea(
                       A2, {CollisionDetection::Segment(A, B),
                            CollisionDetection::Segment(A, D),
                            CollisionDetection::Segment(D, C)}),
                   CollisionDetection::CollisionArea(
                       A2, {CollisionDetection::Segment(B, C)})},
                  {{0, Solid::BB_RIGID_WALL}, {1, Solid::BB_MOVING_WALL}}),

      Config<DIM>({200, 200}, /*iters*/ 30000, /*frames*/ 100,
                  /*reyn*/ 1000.0, /*init_vel*/ {0.1, 0},
                  "out/norms_lid_cavity_openmp_200_1000_01_bgk.bin",
                  "out/data_lid_cavity_openmp_200_1000_01_bgk.bin",
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
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map] = conf;

    types::boundary_mask_t obstacle_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(out_frames);

    Simulation simulation(
        grid_size, obstacle_mask,
        Params<DIM, CollisionType>(reyn, grid_size, init_vel));

    simulation.attachListener(writer);

    MPISolver2D<CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver /*, preconditioner*/, problem);

    simulation.output(out_data.c_str());

    // Confronto con Ghia et al. (1982): lid_velocity = init_vel.dx, la
    // velocita' della parete mobile usata anche per normalizzare i dati
    // tabulati. Richiede reyn in {100, 400, 1000}. Norma scelta qui: L2.
    const auto ghia_re = analysis::ghia_reynolds_from_value(reyn);
    const auto ghia_norm = analysis::NormType::L2;
    const auto ghia = simulation.compute_ghia_error(ghia_re, init_vel.dx, ghia_norm);
    std::cout << "Ghia (" << analysis::to_string(ghia_norm) << ") u: rel="
              << ghia.u_error.relative << " abs=" << ghia.u_error.absolute
              << " | v: rel=" << ghia.v_error.relative
              << " abs=" << ghia.v_error.absolute << std::endl;

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}