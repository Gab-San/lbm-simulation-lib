#include "lbm-2-lbm/types/collision_operators.hpp"
#include "lbm-2-lbm/types/common.hpp"
#include "lbm-2-lbm/core/velocity_sets.hpp"
#include "lbm-2-lbm/solver/solver_2d.hpp"
#include "lbm-2-lbm/lbm_simulation.hpp"
#include "lbm-2-lbm/problems/lid_cavity.hpp"

#include "defines.h"

int main() {
    using namespace lbm;

    LBMSimulation<
    DIM,
    D2Q9,
    MPISolver2D<types::BGK2D>,
    LidCavityProblem<2>,
    CollisionOperatorType::BGK
    > simulation(
	NUM_ITERS, NUM_FRAMES, /* num_skips = num_steps / num_frames */
	INITIAL_VELOCITY,
	types::DimPoint<DIM>(GRID_SIZE_X, GRID_SIZE_Y),
	Params<DIM, CollisionOperatorType::BGK>(
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
    simulation.solve();
    simulation.output("out/data_129_100_01.txt");

    return 0;
}
