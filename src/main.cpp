#include "LBM.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdexcept>

// Max number of steps
#define NSTEPS 5000

// TODO: lid velocity, stream, collide, apply no slip, compute Rho and Velocity, equilibrium,
// TODO: save result + plotting + video/imaging

int main(int argc, char* argv[])
{
    // Command argument count
    if(argc < 5)
    {
        printf("Usage: %s grid_dim_X grid_dim_Y Reyn_Num U_lid\n", argv[0]);
        return -1;
    }

    /// Defining Nx (grid_cells_x) Ny (grid_cells_y)
    /// Re (Reynold number) u_lid (Lid initial velocity)
    int Nx, Ny;
    double Re, u_lid;
    try {
        Nx = std::stoi(argv[1]);    //StringToInteger
        Ny = std::stoi(argv[2]);
        Re = std::stod(argv[3]); //StringToDouble
        u_lid = std::stod(argv[4]);
    } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR WHILE PARSING ARGUMENTS FROM COMMAND LINE] Invalid Argument: " << ia.what() << '\n' << std::endl;
    } catch (const std::out_of_range& oor) {
        std::cout << "[ERROR WHILE PARSING ARGUMENTS FROM COMMAND LINE] Argument Out of Range: " << oor.what() << "\nThis value cannot be represented" << std::endl;
    }

    // Print parameters to stdout
    printf("Simulation parameters:\n");
    printf(" Grid dimensions: %d x %d\n", Nx, Ny);
    printf(" Reynolds number: %f\n", Re);
    printf(" Lid velocity: %f\n", u_lid);

    // Instantiate LBM data structure
    LBM lbm(Nx, Ny, Re, u_lid);

    // definisco le dimensioni da allocare
    size_t mem_size_scalar = Nx * Ny * sizeof(double);        // for rho, ux, uy
    size_t mem_size_ndir   = Nx * Ny * lbm.ndir * sizeof(double);    // for f1, f2

    // allocate memory
    double *f1 = (double*) malloc(mem_size_ndir);
    double *f2 = (double*) malloc(mem_size_ndir);

    // allocate memory for the density field
    double *rho = (double*) malloc(mem_size_scalar); 

    // allocate memory for the velocity fields
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

    // Write grid dimensions to file
    file << Nx << "\n" << Ny << "\n";
    
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {
        // stream from f1 storing to f2
        lbm.stream(f1,f2);
        
        // calculate post-streaming density and velocity
        lbm.compute_rho_u(f2,rho,ux,uy);

        // perform collision on f2
        lbm.collide(f2,rho,ux,uy);

        // apply boundary conditions
        lbm.apply_boundary_conditions(ux,uy,rho,f1,f2);

        // swap pointers
        double *temp = f1;
        f1 = f2;
        f2 = temp;
        
        if(n % 50 != 0) { continue; } // write every 50 time steps

        // calcolo le norme delle velocità e le salvo nel file
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double vx = ux[Nx * j + i];
                double vy = uy[Nx * j + i];
                double v  = sqrt(vx * vx + vy * vy);
                file << v << "\n";
            }
        }
    }

    file.close();
    printf("\nFinished writing to file.\n");

    // deallocate memory
    free(f1); free(f2);
    free(rho); free(ux); free(uy);
    return 0;
}
