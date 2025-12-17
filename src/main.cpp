#include "LBM.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <omp.h>

// Max number of steps
#define NSTEPS 20000

// Skip step
# define SKIP_STEP 100

int main(int argc, char* argv[])
{
    // Command argument count
    if(argc < 5)
    {
        printf("Usage: %s grid_cells_num_x grid_cells_num_y reyn_num lid_velocity\n", argv[0]);
        return -1;
    }

    /// Defining Nx (grid_cells_x) Ny (grid_cells_y)
    int Nx, Ny;
    /// Re (Reynold number) u_lid (Lid initial velocity)
    double Re, u_lid;

    try {
        Nx = std::stoi(argv[1]); //StringToInteger
        Ny = std::stoi(argv[2]);
        Re = std::stod(argv[3]); //StringToDouble
        u_lid = std::stod(argv[4]);
    } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR WHILE PARSING ARGUMENTS FROM COMMAND LINE] Invalid Argument: " 
	    << ia.what() << '\n' << std::endl;
    } catch (const std::out_of_range& oor) {
        std::cout << "[ERROR WHILE PARSING ARGUMENTS FROM COMMAND LINE] Argument Out of Range: " 
	    << oor.what() << "\nThis value cannot be represented" << std::endl;
    }

    // Print parameters to stdout
    printf("Simulation parameters:\n");
    printf(" Grid dimensions: %d x %d\n", Nx, Ny);
    printf(" Reynolds number: %f\n", Re);
    printf(" Lid velocity: %f\n", u_lid);

    #pragma omp parallel
    {
        int n_threads = omp_get_num_threads();
        #pragma omp single
        {
            printf(" Running with %d threads\n", n_threads);
        }
    } 

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
    std::ofstream file("../vel_norms.txt");

    if (!file.is_open()) {
        std::cerr << "Error with opening the file 'vel_norms.txt'.\n";
        return -1;
    }

    // Apertura file per il testing
    std::ofstream fileTest("../vel_y_testing.txt");
    if (!fileTest.is_open()) {
        std::cerr << "Error with opening the file 'vel_y_testing.txt'.\n";
        return -1;
    }

    // Write grid dimensions to file
    file << Nx << "\n" << Ny << "\n";
    
    // main simulation loop; take NSTEPS time steps
    // also compute the time taken for the simulation

    auto startTime = std::chrono::high_resolution_clock::now();

    for(unsigned int n = 0; n < NSTEPS; ++n)
    {
        // stream from f1 storing to f2
        //lbm.stream(f1,f2);

        // apply boundary conditions
        //lbm.apply_boundary_conditions(f2);
        
        // calculate post-streaming density and velocity
        //lbm.compute_rho_u(f2,rho,ux,uy);

        // perform collision on f2
        //lbm.collide(f2,rho,ux,uy);

        lbm.update_stream_collide(f1, f2, rho, ux, uy);

        // swap pointers
        std::swap(f1,f2);
        
        if(n % SKIP_STEP != 0) { continue; } // write every 50 time steps

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
        if (n +  SKIP_STEP >= NSTEPS ) {       
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
    return 0;
}
