#include "LBM.hpp"
#include "config_file_parser/ConfigParser.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <omp.h>

#define WEAK_SCALING false

/**
 * \mainpage Lattice Boltzmann Method Implementation
 *
 * \section Introduction
 *
 * This is an implementation of the lattice Boltzmann method.
 *
 * \section Collision Operators
 *
 * This implementations supports only BGK collision operator.
 *
 * \subsection TRT
 *
 * Using the TRT operator the code would run correctly 
 * but doesn't correctly implement the TRT operator.
 */
int main(int argc, char* argv[])
{
    // Command argument count
    if(argc < 2)
    {
        printf("Usage: %s config_file\n", argv[0]);
        return -1;
    }

    using namespace lbm_lbm;

    std::vector<Job> jobs = parse_file(std::string(argv[1])); 

    std::cout << "Number of simulations: " << jobs.size() << std::endl;

    #if WEAK_SCALING
    std::ofstream weak_scaling_data("../weak_scaling_data.txt", std::ios::binary);
    #endif

    unsigned int Nx, Ny;
    double Re, u_lid;
    unsigned int nskips , nsteps;
    std::string filename_norms, filename_bench;
    for(Job& j : jobs) {

	Nx = j.grid_num_cells_x;
	Ny = j.grid_num_cells_y;
	Re = j.reyn_num;
	u_lid = j.u_lid;
	nsteps = j.num_steps;
	nskips = j.num_frames == 0 ? j.num_steps : j.num_steps / j.num_frames;
	filename_norms = j.vel_norm_out;
	filename_bench = j.vel_bench_out;
	
	// Print parameters to stdout
	printf("Simulation parameters:\n");
	printf(" Grid dimensions: %d x %d\n", Nx, Ny); printf(" Reynolds number: %f\n", Re);
	printf(" Lid velocity: %f\n", u_lid);
	printf(" Number of steps: %d\n Number of frames: %d\n", nsteps, j.num_frames);
	size_t mem_size_ndir   = Nx * Ny * LBM::ndir * sizeof(double); // for f1, f2

	// Instantiate LBM data structure
	LBM lbm(Nx, Ny, Re, u_lid, j.coll_op);

	// initialize lid-driven cavity flow
	lbm.init_lid_driven_cavity();

	// initialise f1 as equilibrium for rho, ux, uy
	lbm.init_equilibrium();

	// Apertura file output
	std::cout << "Opening " << filename_norms << std::endl;
	std::ofstream norms_file(filename_norms, std::ios::binary);

	if (!norms_file.is_open()) {
	    std::cerr << "Error with opening the file " << filename_norms << std::endl;
	    return -1;
	}

	// Apertura file per il testing
	std::cout << "Opening " << j.vel_bench_out << std::endl;
	std::ofstream bench_file(j.vel_bench_out, std::ios::binary);

	if (!bench_file.is_open()) {
	    std::cerr << "Error with opening the file " << j.vel_bench_out << std::endl;
	    return -1;
	}

	// Write grid dimensions to file
	norms_file.write(reinterpret_cast<char *>(&Nx), sizeof(int));
	norms_file.write(reinterpret_cast<char *>(&Ny), sizeof(int));

	// main simulation loop; take NSTEPS time steps
	// also compute the time taken for the simulation

	#if WEAK_SCALING
	auto startTime = std::chrono::high_resolution_clock::now();
	#endif

	for(unsigned int n = 0; n < nsteps; ++n)
	{

	    bool save = n % nskips == 0;
	    lbm.update_stream_collide(save);

	    // swap pointers
	    std::swap(lbm.f1,lbm.f2);

	    if(save)  lbm.write_norms(norms_file);

	    // at the last step 
	    if (n + 1 >= nsteps) lbm.write_bench_data(bench_file);
	}

	#if WEAK_SCALING
	auto endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedSeconds = endTime - startTime;

	printf("Simulation completed in %f seconds.\n", elapsedSeconds.count());
	#endif	

	norms_file.close();
	bench_file.close();
	printf("\nFinished writing to file.\n");    

	#if WEAK_SCALING
	std::string output = std::to_string(Nx) + " " + std::to_string(elapsedSeconds.count());
	weak_scaling_data.write(output.c_str(), sizeof(output.c_str()));
	weak_scaling_data.close();
	#endif

    }

    return 0;
}
