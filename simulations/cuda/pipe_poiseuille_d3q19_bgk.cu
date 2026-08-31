// LBM SIM LIB
#include "lbm-sim/analysis/exact-solution.hpp"
#include "lbm-sim/boundaries/utils.hpp"
#include "lbm-sim/collision-operators/collision-params.hpp"
#include "lbm-sim/core/velocity-sets.hpp"
#include "lbm-sim/data/async-binary-writer.hpp"
#include "lbm-sim/functions.hpp"
#include "lbm-sim/lbm-simulation.hpp"
#include "lbm-sim/logging.hpp"
#include "lbm-sim/solver/cuda-solver.cuh"
#include "lbm/config/config-parser.hpp"

// C++ STD LIB
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static constexpr lbm::types::dim_t DIM = 3;
constexpr auto COLLISION = lbm::CollisionModel::BGK;
// NOTE: no ExecutionBackend constant here. This translation unit hard-codes
// CUDASolver, so an OPEN_MP tag would only be a lie the compiler cannot catch.

// Solid layers left between the outermost fluid node and the domain box. Two
// layers guarantee that every direction leaving a near-wall fluid node lands
// on a solid node, so the four lateral domain BCs are genuinely unreachable.
static constexpr int WALL_MARGIN = 2;

// Thresholds used only for the start-up sanity report.
static constexpr double MAX_MACH = 0.1;          // compressibility error ~ Ma^2
static constexpr double MAX_DENSITY_JUMP = 0.05; // rel. rho drop inlet->outlet
static constexpr double TAU_MIN = 0.55;
static constexpr double TAU_MAX = 1.2;
// Halfway bounce-back sits exactly on the wall when the magic parameter
// Lambda = (1/w+ - 1/2)(1/w- - 1/2) equals 3/16. BGK has a single rate, so
// Lambda = (tau - 1/2)^2 and the wall is exact only at tau = 0.9330.
static constexpr double LAMBDA_OPT = 3.0 / 16.0;
// lambda_1^2, lambda_1 = 2.4048 being the first zero of J0: sets the slowest
// mode of the start-up transient, u -> u_ss like exp(-lambda_1^2 nu t / R^2).
static constexpr double J0_FIRST_ZERO_SQ = 5.7832;

int main(int argc, char **argv) {
  using namespace lbm;

  using types::Coordinate;
  using types::DimPoint;

  // --- 1. CONFIGURATION --------------------------------------------------
  if (argc < 2) {
    config::print_usage(argv[0]);
    return 1;
  }

  std::vector<config::SimulationConfig<DIM>> configs;
  try {
    configs = config::parse_config<DIM>(argv[1]);
  } catch (const config::ConfigError &err) {
    std::cerr << "Configuration error: " << err.what() << "\n";
    return 1;
  }
  // --- 2. INSTANTIATE LOGGER ---------------------------------------------
  logging::setup();
  logging::Logger *main_logger = logging::create_or_get_logger("main");

  for (const auto &cfg : configs) {
    const DimPoint<DIM> grid_size(cfg.grid_size);
    const utils::Vector<double, DIM> u0(cfg.u0);

    const int nx = static_cast<int>(grid_size.x);
    const int ny = static_cast<int>(grid_size.y);
    const int nz = static_cast<int>(grid_size.z);

    // --- 2b. GEOMETRY CHECKS ---------------------------------------------
    // Odd ny and nz put the axis exactly on a node column, so the wall is at
    // the same distance on both sides and the centre profile really is a
    // diameter. With an even count the axis falls between two nodes and the
    // parabola comes out asymmetric by half a lattice unit.
    if (ny % 2 == 0 || nz % 2 == 0) {
      std::cerr << "Simulation '" << cfg.name
                << "': ny and nz must be odd for the pipe axis to lie on a "
                   "node column; got "
                << ny << " x " << nz << ". Skipped.\n";
      continue;
    }
    if (ny != nz) {
      std::cerr << "Simulation '" << cfg.name
                << "': ny != nz, the cross-section is not square; the pipe is "
                   "still circular but the margins differ.\n";
    }

    // Pipe axis: the line (y, z) = (cy, cz), parallel to x.
    const int cy = ny / 2;
    const int cz = nz / 2;

    // The radius is derived from the grid instead of being hard-coded, so a
    // config that changes grid_size stays consistent (the old fixed 30 gave a
    // silently wrong geometry for any ny other than 65).
    const int radius = std::min(ny, nz) / 2 - WALL_MARGIN;
    if (radius < 4) {
      std::cerr << "Simulation '" << cfg.name << "': grid too small, derived "
                << "radius = " << radius << ". Skipped.\n";
      continue;
    }

    // With halfway bounce-back the wall sits midway between the last fluid
    // node and the first solid node.
    // CHECK ONCE against your CylindricalShell convention: this assumes the
    // inner radius is exclusive, i.e. nodes at exactly r = radius are still
    // fluid. If they are solid, R_eff is radius - 0.5 and every error below
    // is biased by one lattice unit on the diameter.
    const double r_eff = radius + 0.5;
    const double d_eff = 2.0 * r_eff;

    // --- 3. TRANSPORT PARAMETERS ------------------------------------------
    // CollisionParams derives nu from the characteristic length num_cells.y.
    // For a pipe the characteristic length is the diameter, not the box side,
    // so the Reynolds number handed to it is rescaled: passing Re * ny / D
    // yields nu = u * D / Re, i.e. the Reynolds number actually simulated is
    // the one written in the config. (The clean fix is an overload of
    // CollisionParams taking an explicit characteristic length.)
    const double re_target = cfg.reynolds;
    const double re_scaled = re_target * static_cast<double>(ny) / d_eff;
    const CollisionParams<DIM, COLLISION> params(re_scaled, grid_size, u0);

    const double u_max = params.init_vel.dx;
    const double tau = numbers::invcs_2 * params.nu + 0.5;
    const double mach = std::sqrt(numbers::invcs_2) * u_max;
    const double lambda = (tau - 0.5) * (tau - 0.5);

    // The pressure drop that sustains the flow: as in the 2D channel it is not
    // a free parameter, it is fixed by inverting the Hagen-Poiseuille solution
    // to give u_max = init_vel.dx on the axis. In lattice units p = rho*cs^2,
    // so the *density* jump is invcs_2 times the pressure one; the two
    // arguments below are therefore densities, hence the names.
    //
    // The driving length is the period of the pressure-periodic pair, nx.
    // If your implementation puts inlet and outlet on nodes 0 and nx-1 with no
    // wrap-around cell, use nx - 1 here; the difference is O(1/nx).
    const double length = static_cast<double>(nx);
    const double delta_rho =
        numbers::invcs_2 * 4.0 * params.nu * u_max * length / (r_eff * r_eff);
    const double rho_out = 1.0;
    const double rho_in = rho_out + delta_rho;

    // Slowest mode of the start-up transient. Reaching 0.1% of steady state
    // takes about 7 of these.
    const double t_transient = (r_eff * r_eff) / (J0_FIRST_ZERO_SQ * params.nu);

    LBM_LOG_INFO(
        main_logger,
        "Simulation '{}':\n\tGrid dimensions: {}\n\tPipe axis: (y,z) = "
        "({},{}), radius {} (R_eff {}, D_eff {})\n\tReynolds number: {} (on "
        "D_eff; passed to CollisionParams as {})\n\tReference velocity: "
        "{}\n\tnu: {}, tau: {}, Mach: {}\n\tLambda: {} (optimum {})\n\tDensity "
        "jump: {} ({:.3f}%)\n\tTransient time scale: {:.0f} steps\n\tNumber of "
        "Iterations: {}\n\tNumber of frames: {}\n\tFrames output: "
        "{}\n\tProfile output: {}",
        cfg.name, grid_size, cy, cz, radius, r_eff, d_eff, re_target, re_scaled,
        u0, params.nu, tau, mach, lambda, LAMBDA_OPT, delta_rho,
        100.0 * delta_rho / rho_out, t_transient, cfg.niters, cfg.nframes,
        cfg.frames_out, cfg.profile_out);

    if (tau < TAU_MIN || tau > TAU_MAX) {
      LBM_LOG_NOTICE(main_logger,
                     "'{}': tau = {} outside [{}, {}]: expect poor accuracy or "
                     "instability. Raise the resolution or lower Re.",
                     cfg.name, tau, TAU_MIN, TAU_MAX);
    }
    if (mach > MAX_MACH) {
      LBM_LOG_NOTICE(main_logger,
                     "'{}': Mach = {} > {}: compressibility error is no longer "
                     "negligible. Lower u0.",
                     cfg.name, mach, MAX_MACH);
    }
    if (delta_rho / rho_out > MAX_DENSITY_JUMP) {
      LBM_LOG_NOTICE(main_logger,
                     "'{}': density jump {} of rho_out: the flow is not "
                     "incompressible along x. Shorten the pipe or lower u0.",
                     cfg.name, delta_rho / rho_out);
    }
    if (static_cast<double>(cfg.niters) < 7.0 * t_transient) {
      LBM_LOG_NOTICE(main_logger,
                     "'{}': niters = {} is below the ~{:.0f} steps needed to "
                     "reach steady state: the L2 error will be dominated by "
                     "the transient, not by the discretisation.",
                     cfg.name, cfg.niters, 7.0 * t_transient);
    }

    // --- 4. CREATE OBSTACLES ----------------------------------------------
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
    // included. Both radii are int so the shell's template argument deduces
    // from a single type.
    const CollisionDetection::CollisionArea<DIM> pipe_wall(
        Coordinate<DIM>(0, 0, 0),
        {CollisionDetection::CylindricalShell<DIM>(
            Coordinate<DIM>(0, cy, cz), Coordinate<DIM>(nx - 1, cy, cz), radius,
            ny + nz)});

    // id 0 = the pipe wall: rigid bounce-back, stationary.
    const std::vector<Solid::ObstacleData<DIM>> obstacle_data{
        {Solid::BB_RIGID_WALL, {0.0, 0.0, 0.0}}};

    // --- 5. CREATE MASK ----------------------------------------------------
    types::solid_mask_t solid_mask =
        Solid::compute_solid_mask<DIM>({pipe_wall}, grid_size);

    // --- 6. RUN SIMULATION -------------------------------------------------
    // frames_out is the DIRECTORY; the file basename comes from the config
    // name, so different runs in the same directory do not overwrite each
    // other.
    std::shared_ptr<AsyncBinaryWriter> writer =
        std::make_shared<AsyncBinaryWriter>(cfg.frames_out);

    LBMSimulation<DIM, D3Q19, COLLISION> simulation(
        grid_size, std::move(solid_mask), obstacle_data, dbc, params, rho_in,
        rho_out);
    simulation.attachListener(writer);

    CUDASolver<DIM, D3Q19, COLLISION> solver(cfg.niters, cfg.nframes);
    solver.attachListener(writer);

    simulation.solve(solver);

    // --- 7. OUTPUT ---------------------------------------------------------
    // ux along z on the centre column (x = nx/2, y = ny/2): that is a diameter
    // of the pipe, i.e. the section on which the parabola is read off.
    simulation.output(cfg.profile_out.c_str(),
                      functional::extract_dx_profile_along_z_center);

    // --- 8. ERROR COMPUTATION ----------------------------------------------
    const auto exact_solution =
        analysis::HagenPoiseuilleSolution3D(r_eff, u_max, cy, cz);
    const double err_l2 =
        simulation.compute_error(analysis::NormType::L2, exact_solution);

    LBM_LOG_NOTICE(main_logger, "{} error: {}",
                   analysis::to_string(analysis::NormType::L2), err_l2);

    simulation.detachListener(writer);
    solver.detachListener(writer);
  }
  return 0;
}
