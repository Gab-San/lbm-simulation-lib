// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/solver-2d.hpp"

#include "lbm-sim/collision-detection/collision-area.hpp"

#include "lbm-sim/analysis/exact-solution.hpp"

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

  /// Quale soluzione analitica usare per compute_error() a fine run.
  /// Fissato qui una volta sola: niente H/Umax duplicati altrove.
  const lbm::analysis::FlowType flow_type = lbm::analysis::FlowType::Couette;

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
  using utils::Vector;

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(128, 128);
  const Coordinate<2> D(128, 0);

  std::vector<Config<2>> configs{
      Config<2>(
          {129, 129}, /*iters*/ 100000, /*frames*/ 1000, /*reyn*/ 100.0,
          /*init_vel*/ {0.1, 0}, "out/norms_couette_129_100_01.bin",
          "out/data_couette_129_100_01.bin",
          {
              CollisionDetection::CollisionArea(
                  A, {CollisionDetection::Segment(A, D)}), // bottom (y=0)
              CollisionDetection::CollisionArea(
                  A, {CollisionDetection::Segment(B, C)}), // top (y=128)
              CollisionDetection::CollisionArea(
                  A, {CollisionDetection::Segment(A + Vector<int, DIM>(0, 1),
                                                  B - Vector<int, DIM>(0, 1)),
                      CollisionDetection::Segment(C - Vector<int, DIM>(0, 1),
                                                  D + Vector<int, DIM>(0, 1))}),
          },
          {{0, Solid::BB_RIGID_WALL},
           {1, Solid::BB_MOVING_WALL},
           {2, Solid::PERIODIC}}),
  };

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;
  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map, flow_type] = conf;
    types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(conf.out_frames);

    Simulation simulation(
        grid_size, boundary_mask,
        Params<DIM, CollisionType>(reyn, grid_size, init_vel));

    simulation.attachListener(writer);

    MPISolver2D<CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver, problem);
    simulation.output(out_data.c_str());

    // H = altezza canale (parete inferiore a y=0, superiore a y=grid_size.y-1);
    // Umax = velocita' di riferimento (parete mobile per Couette).
    // Stessi valori gia' usati per costruire la simulazione: nessuna
    // duplicazione, flow_type sceglie la Function<2> corretta.
    const double H = static_cast<double>(grid_size.y - 1);
    const auto exact_solution =
        analysis::make_exact_solution(flow_type, H, init_vel.dx);
    const double err_l2 =
        simulation.compute_error(analysis::NormType::L2, *exact_solution);
    std::cout << "Errore L2 vs soluzione analitica: " << err_l2 << std::endl;

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}