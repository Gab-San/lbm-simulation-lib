// Pipe flow: 3D Hagen-Poiseuille in a duct of circular cross-section,
// D3Q19 velocity set, BGK operator, OpenMP backend.
//
// This is the 3D sibling of poiseuille_flow_2d_bgk.cpp: same engine
// (pressure imposed at inlet and outlet, bounce-back on the wall), but the
// plane channel between two walls becomes a cylindrical duct. The expected
// profile changes accordingly -- a paraboloid of revolution instead of a
// parabola in y -- and so does the pressure drop that sustains it:
//
//   channel (2D):  u_max = dp * H^2 / (8 * mu * L)
//   pipe    (3D):  u_max = dp * R^2 / (4 * mu * L)      <- Hagen-Poiseuille
//
// The pipe wall is not a domain face but an immersed obstacle: the domain
// stays a box and the CylindricalShell declares everything beyond the pipe
// radius solid, carving the duct out of the box. The box corners stay solid
// and take no part in the flow.

// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/collision-detection/collision-area.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/config/simulation-config.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/vtk-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/openmp-solver.hpp"

// C++ STD LIB
#include <memory>
#include <string>
#include <vector>

static constexpr lbm::types::dim_t DIM = 3;
constexpr auto COLLISION = lbm::CollisionModel::BGK;
constexpr auto BACKEND = lbm::ExecutionBackend::OPEN_MP;

int main(int argc, char **argv) {
  using namespace lbm;

  using types::Coordinate;
  using types::DimPoint;
  using utils::Vector;

  // --- 1. CONFIGURATION --------------------------------------------------
  config::SimulationConfig<DIM> cfg;

  cfg.name = "pipe_flow_3d_bgk";

  cfg.backend = BACKEND;
  cfg.collision = COLLISION;

  // Pipe along x. The cross-section is square and the duct fits inside it
  // with two cells of margin per side: the wall therefore falls *inside*
  // the domain, where bounce-back can represent it, and not on the faces.
  cfg.grid_size = {128, 65, 65};

  cfg.u0 = {0.1, 0, 0};

  cfg.reynolds = 100;

  cfg.niters = 15000;
  cfg.nframes = 300;

  cfg.frames_out = "output/pipe_flow_3d_bgk_frames";
  cfg.profile_out = "output/pipe_flow_3d_bgk_profile.dat";

  // --- 2. INSTANTIATE LOGGER ---------------------------------------------
  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  const DimPoint<DIM> grid_size(cfg.grid_size);

  const int nx = static_cast<int>(grid_size.x);
  const int ny = static_cast<int>(grid_size.y);
  const int nz = static_cast<int>(grid_size.z);

  // Pipe axis: the line (y, z) = (cy, cz), parallel to x.
  const int cy = ny / 2;
  const int cz = nz / 2;
  const unsigned int radius = 30;

  LBM_LOG_INFO(
      main_logger,
      "Simulation '{}':\n\tGrid dimensions: {}\n\tPipe axis: (y,z) = "
      "({},{}), radius {}\n\tReynolds number: {}\n\tReference velocity: "
      "{}\n\tNumber of Iterations: {}\n\tNumber of frames: {}\n\tFrames "
      "output: {}\n\tProfile output: {}",
      cfg.name, grid_size, cy, cz, radius, cfg.reynolds, cfg.u0, cfg.niters,
      cfg.nframes, cfg.frames_out, cfg.profile_out);

  // --- 3. CREATE OBSTACLES -----------------------------------------------
  // Inlet and outlet are pressure-imposed on the two x faces; the other four
  // faces are already buried in the solid of the wall, so the BC they carry
  // is never reached by a fluid node: rigid is the consistent choice anyway.
  Solid::DomainBC<DIM> dbc{};
  dbc.low(0) = Solid::PRESSURE_PERIODIC_INLET;   // x = 0
  dbc.high(0) = Solid::PRESSURE_PERIODIC_OUTLET; // x = nx-1
  dbc.low(1) = Solid::BB_RIGID_WALL;             // y = 0
  dbc.high(1) = Solid::BB_RIGID_WALL;            // y = ny-1
  dbc.low(2) = Solid::BB_RIGID_WALL;             // z = 0
  dbc.high(2) = Solid::BB_RIGID_WALL;            // z = nz-1

  // The pipe wall: a cylindrical shell from x = 0 to x = nx-1 around the
  // axis. The outer radius (ny + nz) is not a real radius, it just means
  // "large enough to leave the domain": the AABB gets clipped to the grid,
  // so the shell in practice covers every node beyond `radius`, box corners
  // included.
  const CollisionDetection::CollisionArea<DIM> pipe_wall(
      Coordinate<DIM>(0, 0, 0),
      {CollisionDetection::CylindricalShell<DIM>(
          Coordinate<DIM>(0, cy, cz), Coordinate<DIM>(nx - 1, cy, cz), radius,
          ny + nz)});

  // id 0 = the pipe wall: rigid bounce-back, stationary.
  const std::vector<Solid::ObstacleData<DIM>> obstacle_data{
      {Solid::BB_RIGID_WALL, {0.0, 0.0, 0.0}}};

  // --- 4. CREATE MASK ----------------------------------------------------
  types::solid_mask_t solid_mask =
      Solid::compute_solid_mask<DIM>({pipe_wall}, grid_size);

  // --- 5. RUN SIMULATION -------------------------------------------------
  // frames_out is the DIRECTORY; the file basename comes from the config
  // name, so different runs in the same directory do not overwrite each
  // other.
  std::shared_ptr<VtkWriter> writer =
      std::make_shared<VtkWriter>(cfg.frames_out, cfg.name);

  // WARNING: CollisionParams derives nu from the characteristic length
  // num_cells.y, not from the pipe diameter. With ny = 65 and an effective
  // diameter 2R+1 = 61 the Reynolds number actually simulated is
  // u*D/nu ~ 94, not 100: keep the radius close to ny/2 if you want the two
  // to agree.
  const CollisionParams<DIM, COLLISION> params(cfg.reynolds, grid_size, cfg.u0);

  // The pressure drop that sustains the flow: as in the 2D channel it is not
  // a free parameter, it is fixed by inverting the Hagen-Poiseuille solution
  // to give u_max = init_vel.dx on the axis. In lattice units p = rho*cs^2,
  // so the *density* jump is invcs_2 times the pressure one.
  //
  // R is the effective wall radius: with halfway bounce-back the wall sits
  // midway between the last fluid node and the first solid node.
  const double R_eff = radius + 0.5;
  const double pout = 1;
  const double pin = pout + numbers::invcs_2 * 4.0 * params.nu *
                                params.init_vel.dx * nx / (R_eff * R_eff);

  LBMSimulation<DIM, D3Q19, COLLISION> simulation(
      grid_size, std::move(solid_mask), obstacle_data, dbc, params, pin, pout);
  simulation.attachListener(writer);

  OpenMPSolver<DIM, D3Q19, COLLISION> solver(cfg.niters, cfg.nframes);
  solver.attachListener(writer);

  simulation.solve(solver);

  // --- 6. OUTPUT ---------------------------------------------------------
  // ux along z on the centre column (x = nx/2, y = ny/2): that is a diameter
  // of the pipe, i.e. the section on which the parabola is read off.
  simulation.output(cfg.profile_out.c_str(),
                    functional::extract_dx_profile_along_z_center);

  // --- 7. ERROR COMPUTATION ----------------------------------------------
  const auto exact_solution =
      analysis::HagenPoiseuilleSolution3D(R_eff, cfg.u0.dx, cy, cz);
  const double err_l2 =
      simulation.compute_error(analysis::NormType::L2, exact_solution);

  LBM_LOG_NOTICE(main_logger, "{} error: {}",
                 analysis::to_string(analysis::NormType::L2), err_l2);

  simulation.detachListener(writer);
  solver.detachListener(writer);

  return 0;
}
