#include "LBM.hpp"
#include "config_file_parser/ConfigParser.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <omp.h>

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
    
    std::cout << "Number of jobs: " << jobs.size() << std::endl;

    unsigned int Nx, Ny;
    double Re, u_lid;
    unsigned int nskips , nsteps;
    std::string filename;
    for(Job& j : jobs) {

	Nx = j.grid_num_cells_x;
	Ny = j.grid_num_cells_y;
	Re = j.reyn_num;
	u_lid = j.u_lid;
	nsteps = j.num_steps;
	nskips = j.num_steps / j.num_frames;
	filename = j.vel_norm_out;

	// Print parameters to stdout
	printf("Simulation parameters:\n");
	printf(" Grid dimensions: %d x %d\n", Nx, Ny);
	printf(" Reynolds number: %f\n", Re);
	printf(" Lid velocity: %f\n", u_lid);

	// Instantiate LBM data structure
	LBM lbm(Nx, Ny, Re, u_lid);

	// Allocation sizes
	// for scalar valued vectors
	size_t mem_size_scalar = Nx * Ny * sizeof(double);            // for rho, ux, uy
	// for vectors that take into account the neighbouring cells
	size_t mem_size_ndir   = Nx * Ny * lbm.ndir * sizeof(double); // for f1, f2

	// allocate memory
	double *f1 = (double*) malloc(mem_size_ndir);
	double *f2 = (double*) malloc(mem_size_ndir);

	// allocate memory for the density field
	double *rho = (double*) malloc(mem_size_scalar); 

	// allocate memory for the velocity fields
	// considering that the problem is the lid cavity problem
	// an initial optimization can be done by setting both 
	// velocity vectors to 0
	double *ux = (double*) malloc(mem_size_scalar);
	double *uy = (double*) malloc(mem_size_scalar);

	// compute Taylor-Green flow at t=0
	// to initialise rho, ux, uy fields.
	//lbm.taylor_green(0,rho,ux,uy);

	// initialize lid-driven cavity flow
	lbm.init_lid_driven_cavity(ux, uy, rho);

	// initialise f1 as equilibrium for rho, ux, uy
	lbm.init_equilibrium(f1,rho,ux,uy);

	// Apertura file output
	std::cout << "Opening " << filename << std::endl;
	std::ofstream file(filename);

	if (!file.is_open()) {
	    std::cerr << "Error with opening the file " << filename << std::endl;
	    return -1;
	}

	// Apertura file per il testing
	std::cout << "Opening " << j.vel_bench_out << std::endl;
	std::ofstream fileTest(j.vel_bench_out);

	if (!fileTest.is_open()) {
	    std::cerr << "Error with opening the file " << j.vel_bench_out << std::endl;
	    return -1;
	}

	// Write grid dimensions to file
	file << Nx << "\n" << Ny << "\n";

	// main simulation loop; take NSTEPS time steps
	// also compute the time taken for the simulation

	auto startTime = std::chrono::high_resolution_clock::now();

	for(unsigned int n = 0; n < nsteps; ++n)
	{

	    lbm.update_stream_collide(f1, f2, rho, ux, uy);

	    // swap pointers
	    std::swap(f1,f2);

	    if(n % nskips != 0) { continue; } // write every 50 time steps

	    // saving euclidian norms of the vectors at each step 
	    for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
		    double vx = ux[Nx * j + i];
		    double vy = uy[Nx * j + i];
		    double v  = sqrt(vx * vx + vy * vy);
		    file << v << "\n";
		}
	    }

	    // at the last step 
	    if (n +  nskips >= nsteps ) {       
		// Results for v-Velocity along Horizontal Line through Geometric Center of Cavity
		int j_center = Ny / 2;
		for (int i = 0; i < Nx; ++i) {
		    double v_center = uy[Nx * j_center + i];
		    fileTest << v_center << "\n";
		}
		
	    }
	}

	auto endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedSeconds = endTime - startTime;

	printf("Simulation completed in %f seconds.\n", elapsedSeconds.count());

	file.close();
	fileTest.close();
	printf("\nFinished writing to file.\n");    

	// deallocate memory
	free(f1); free(f2);
	free(rho); free(ux); free(uy);
    }

    return 0;
}
