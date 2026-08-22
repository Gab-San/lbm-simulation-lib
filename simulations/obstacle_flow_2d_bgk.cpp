// LBM SIM LIB
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

#include "lbm-sim/boundaries.hpp"
#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/omp-solver.hpp"

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

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(128, 128);
  const Coordinate<2> D(128, 0);

  std::vector<Config<2>> configs{
      Config<2>(
          {129, 129}, /*iters*/ 100000, /*frames*/ 200, /*reyn*/ 100.0,
          /*init_vel*/ {0.1, 0}, "out/norms_obstacle_129_100_01_bgk.bin",
          "out/data_obstacle_129_100_01_bgk.bin",
          {
              CollisionDetection::CollisionArea(
                  A, {CollisionDetection::Segment(A, D),   // bottom (y=0)
                      CollisionDetection::Segment(B, C)}), // top (y=128)
              CollisionDetection::CollisionArea(           // LEFT WALL
                  A, {CollisionDetection::Segment(A + Vector<int, DIM>(0, 1),
                                                  B - Vector<int, DIM>(0, 1))}),
              CollisionDetection::CollisionArea( // RIGHT WALL
                  A, {CollisionDetection::Segment(C - Vector<int, DIM>(0, 1),
                                                  D + Vector<int, DIM>(0, 1))}),

              // Nuovo ostacolo dentro il dominio
              // CollisionDetection::CollisionArea(
              //     Coordinate<2>(0, 0),         // posizione base (l'offset per le coord del cerchio)
              //     {CollisionDetection::Circle<DIM>(
              //         Coordinate<2>(64, 64),  // centro relativo alla posizione base
              //         16)}                    // raggio in celle
              // ),

              CollisionDetection::CollisionArea(
                    Coordinate<2>(0, 0),
                    {CollisionDetection::Airfoil<DIM>(
                        Coordinate<2>(40, 64),  // leading edge del profilo
                        50.0,                   // corda in celle
                        0.12,                   // spessore         0.12 -> NACA XX12
                        0.02,                   // camber massimo   0.02 -> NACA 2XXX
                        0.4,                    // posizione camber 0.40 -> NACA X4XX
                        15.0)}                   // angolo di attacco
              ),
          },
          {{0, Solid::BB_RIGID_WALL},   // fixed top and bottom wall
           {1, Solid::PRESSURE_PERIODIC_INLET},
           {2, Solid::PRESSURE_PERIODIC_OUTLET}, // right and left
           {3, Solid::BB_OBSTACLE}}), // <= indice 3 = il nuovo oggetto
  };

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  const LidCavity2D problem;

  for (const auto &conf : configs) {
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obst_type_map] = conf;
    types::boundary_mask_t boundary_mask =
        Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(conf.out_frames);

    Params<DIM, CollisionType> params(reyn, grid_size, init_vel);
    const double pout = 1;
    const double pin =
        pout + (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) *
                   8 * params.nu * params.init_vel.dx;
    Simulation simulation(grid_size, boundary_mask, params, pin, pout);

    simulation.attachListener(writer);

    MPISolver2D<CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver, problem);
    simulation.output(out_data.c_str());

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
