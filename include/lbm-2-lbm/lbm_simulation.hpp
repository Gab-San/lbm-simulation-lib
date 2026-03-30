#pragma once

#include "lbm-2-lbm/core/collision_operators.hpp"
#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/types/common.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <filesystem>

namespace lbm {

template<
int dim,
class VelocitySet,
class Solver,
class Problem,
enum CollisionOperatorType cs_type = CollisionOperatorType::BGK
>
class LBMSimulation {
private:
	static_assert(
	dim == VelocitySet::dim,
	"LBMSimulation: template parameter 'dim' must match VelocitySet::dim"
	);
	static_assert(
	cs_type == Solver::CS_TYPE,
	"LBMSimulation: collision operator type must match Solver::CS_TYPE"
	);

    // ----- TEMPLATE VARS ------
    Grid<dim> grid;

    // ----- PASSED VARS -------
    const int niters, nskips;
    const double init_vel;
    const Params<dim, cs_type> params;

public:

    LBMSimulation(
	const int _num_iters,
	const int _num_frames,
	const double _init_vel,
	const types::DimPoint<dim> _grid_dim,
	const Params<dim, cs_type> _params
    ) :
	niters(_num_iters),
	nskips(_num_iters / _num_frames), 
	init_vel(_init_vel),
	params(_params),
	grid(_grid_dim)
    {
	if (_num_iters < _num_frames)
	    throw std::invalid_argument(
		"the number of iterations must be higher than the number of frames!"
	    );
    };

    void setup(
	#if listeners
	FileWriter fw
	#endif
    ) {};

    void solve() {
	std::cout << "Initializing Simulation." << std::endl;

	Problem problem(grid, init_vel);
	Solver solver(grid, init_vel, params);
	std::vector<double> f1(grid.getArea() * VelocitySet::ndir, 0.0);
	std::vector<double> f2(grid.getArea() * VelocitySet::ndir, 0.0);

	problem.init();
	std::cout << "Problem Initialized." << std::endl;
	solver.init_equilibrium(f1);
	std::cout << "Equilibrium Initialized." << std::endl;

	#if listeners
	fw.write(GridSize);
	#endif

	for( unsigned int iter = 0; iter < niters; iter++ ) {
	    bool save = iter % nskips == 0;
	    solver.stream(f1, f2);
	    problem.apply_boundary_conditions(f2);
	    solver.update_collide(f1, f2, save);
		std::swap(f1, f2);
		
	    #if listeners
	    if ( save ) fw.write(WriteNorms);
	    #endif
	}

	std::cout << "Finished Simulation." << std::endl;
    };

    #if listeners
    void notifyAll(Event e) {
	for(Listener listener : simulationListeners) {
	    listener.listen(e);
	}
    }
    #endif

    void output(const char* filepath) {
	using namespace std::filesystem;

	path outpath = filepath;
	path parent = outpath.parent_path();

	if (!exists(parent))
	    create_directories(parent);

	std::cout << "Opening " << outpath << std::endl;

	std::ofstream fout(outpath, std::ios::binary);

	if (!fout.is_open())
	    std::cerr << "Failed to create file: " <<
		outpath << std::endl;

	std::cout << "Writing..." << std::endl;

	std::vector<double> v_center(grid.Nx);
	int j_center = grid.Ny / 2;
	for (int i = 0; i < grid.Nx; ++i) {
	    v_center[i] = grid.uy[grid.Nx * j_center + i];
	}

	fout.write(
	    reinterpret_cast<const char*>(v_center.data()), 
	    v_center.size() * sizeof(double)
	);

	std::cout << "Finished writing to " << outpath << std::endl;

	fout.close();
    };
};

}
