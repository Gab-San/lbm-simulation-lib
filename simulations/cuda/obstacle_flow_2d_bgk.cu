// LBM SIM LIB
#include "lbm-sim/collision-operators/collision-params.hpp"

#include "lbm-sim/core/velocity-sets.hpp"

#include "lbm-sim/data/async-binary-writer.hpp"

#include "lbm-sim/boundaries/boundary-conditions.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/lbm-simulation.hpp"

#include "lbm-sim/solver/cuda-solver.cuh"

#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/functions.hpp"

#include "lbm-sim/logging.hpp"

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

  /// Reference velocity, used only to compute nu/tau.
  /// It no longer moves any wall in Poiseuille!
  const lbm::utils::Vector<double, 2> init_vel;

  /// Output path for frames
  const std::string out_frames;

  /// Output path for benchmark data
  const std::string out_data;

  /// Only bodies immersed in the fluid: the channel walls live in domain_bc.
  const std::vector<lbm::CollisionDetection::CollisionArea<DIM>> obstacles;

  /// Side table: obstacle id -> {BC type, wall velocity}.
  const std::vector<lbm::Solid::ObstacleData<DIM>> obstacle_data;

  /// BCs of the four domain faces.
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

/// Poiseuille channel: pressure imposed at inlet and outlet, rigid top and
/// bottom walls.
static lbm::Solid::DomainBC<DIM> make_channel_bc() {
  lbm::Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = lbm::Solid::PRESSURE_PERIODIC_INLET;   // x = 0
  dbc.high(0) = lbm::Solid::PRESSURE_PERIODIC_OUTLET; // x = nx-1
  dbc.low(1) = lbm::Solid::BB_RIGID_WALL;             // y = 0
  dbc.high(1) = lbm::Solid::BB_RIGID_WALL;            // y = ny-1
  return dbc;
}

int main() {
  using namespace lbm;
  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  const Coordinate<2> A(0, 0);
  const Coordinate<2> B(0, 128);
  const Coordinate<2> C(639, 128);
  const Coordinate<2> D(639, 0);

  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  std::vector<Config<2>> configs{
      Config<2>(
          {640, 129}, /*iters*/ 1000000, /*frames*/ 200, /*reyn*/ 1000.0,
          /*init_vel*/ {0.05, 0}, "out/norms_obstacle_129_100_01_bgk.bin",
          "out/data_obstacle_129_100_01_bgk.bin",
          {
              // The channel walls are no longer obstacles: they live in
              // make_channel_bc(). Only the immersed body is left here.
              // CollisionDetection::CollisionArea(
              //     Coordinate<2>(
              //         0,
              //         0), // base position (offset of the circle coordinates)
              //     {CollisionDetection::Circle<DIM>(
              //         Coordinate<2>(160,
              //                       64), // centre relative to base position
              //         16)}               // radius in cells
              // ),

              // NACA airfoil profile
              CollisionDetection::CollisionArea(
                  Coordinate<2>(0, 0),
                  {CollisionDetection::Airfoil<DIM>(
                      Coordinate<2>(100, 64), // leading edge of the profile
                      100.0,                  // chord, in cells
                      0.12,                   // thickness    0.12 -> NACA XX12
                      0.02,                   // max camber   0.02 -> NACA 2XXX
                      0.4,                    // camber pos.  0.40 -> NACA X4XX
                      5.0)}                   // angle of attack
                  ),

          },
          // id 0 = the cylinder: rigid wall, stationary.
          {{Solid::BB_RIGID_WALL, {0.0, 0.0}}}, make_channel_bc()),
  };

  constexpr auto CollisionType = CollisionModel::BGK;
  using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

  for (auto &conf : configs) {
    const auto &[grid_size, iters, frames, reyn, init_vel, out_frames, out_data,
                 obstacles, obstacle_data, domain_bc] = conf;

    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>(obstacles, grid_size);

    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(conf.out_frames);

    CollisionParams<DIM, CollisionType> params(reyn, grid_size, init_vel);
    const double pout = 1;
    const double pin =
        pout +
        numbers::invcs_2 *
            (grid_size.x / static_cast<double>(grid_size.y * grid_size.y)) * 8 *
            params.nu * params.init_vel.dx;

    Simulation simulation(grid_size, std::move(solid_mask), obstacle_data,
                          domain_bc, params, pin, pout);
    simulation.attachListener(writer);

    CUDASolver<DIM, D2Q9, CollisionType> solver(iters, frames);
    solver.attachListener(writer);

    simulation.solve(solver);
    simulation.output(out_data.c_str(),
                      functional::extract_dx_profile_along_y_center);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }

  return 0;
}
