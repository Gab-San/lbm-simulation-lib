#include "lbm-2-lbm/lbm_simulation.hpp"

#include "lbm-2-lbm/core/velocity_sets.hpp"

#include "lbm-2-lbm/types/common.hpp"
#include "lbm-2-lbm/types/collision_operators.hpp"

#include "lbm-2-lbm/solver/solver_2d.hpp"

#include "lbm-2-lbm/problems/problem_2d.hpp"
#include "lbm-2-lbm/problems/naca_airfoil_2d.hpp"


#include "defines.h"

#include <memory>

int main() {
    using namespace lbm;

    //const auto problem = std::make_unique<LidCavity2D<ExecutionBackend::MPI>>();

    const auto problem = std::make_unique<NACAirfoil2D<ExecutionBackend::MPI>>(
        "NACA2412.dat",  // path al file .dat
        0.30,            // chord / Nx
        0.25,            // leading edge x / Nx
        0.50,            // centro verticale / Ny
        5.0,             // angle of attack [deg]
        true             // open channel (free-slip top/bottom)
    );

    LBMSimulation<
    DIM,
    D2Q9,
    MPISolver2D<types::BGK2D>,
    types::BGK,
    ExecutionBackend::MPI
    > simulation(
        NUM_ITERS, NUM_FRAMES, /* num_skips = num_steps / num_frames */
        INITIAL_VELOCITY,
        types::DimPoint<DIM>(GRID_SIZE_X, GRID_SIZE_Y),
        Params<DIM, types::BGK>(
            REYNOLDS_NUM, GRID_SIZE_Y
        )
    );

    simulation.setup(
        #if listeners
        SimOutWriter<NORM>(),
        SimOutWriter<BENCH, X>(),
        SimOutWriter<BENCH, Y>(),
        SimOutWriter<BENCH, Z>()
        #endif
    );
    const ProblemBase<DIM, ExecutionBackend::MPI>& selected_problem = *problem;
    simulation.solve(selected_problem);
    simulation.output("out/data_129_100_01.txt");

    return 0;
}
