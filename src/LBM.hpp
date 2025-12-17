#endif
#ifndef LBM_HPP
#define LBM_HPP

#include <stdexcept> // for runtime_error
#include <cstddef>   // for size_t
#include <vector>

/**
 * LBM
 *
 * This class handles the functions through which the lattice Boltzmann
 * method can be executed.
 *
 * It contains the general parameters needed to implement the algorithm
 * and it is initialized with the parameters of the problem.
 */
class LBM {

public:

    /**
     * Default Construct
     *
     * @param grid_x_ number of cells along the X axis
     * @param grid_y_ number of cells along the Y axis
     * @param rey_num_ : Reynold number (viscosity related)
     * @param u_lid : Initial velocity of the lid cavity
     * @param nu: kinematic viscosity
     * @param tau: relaxation parameter
     */
    LBM(int num_cells_x_, int num_cells_y_, double rey_num_, double u_lid_):
        nu(u_lid_ * (num_cells_y_) / rey_num_), 
        tau(0.5 + 3.0 * nu),
        Nx(num_cells_x_), 
        Ny(num_cells_y_), 
        u_lid(u_lid_),
        Re(rey_num_),
        taup_inv(2.0/(6.0*nu+1.0)),
        taum_inv(2.0/(6.0*nu+1.0))
    /*
    two tau relaxation for TRT
    */
        {
        if (tau <= 0.5) {
            throw std::runtime_error("LBM error: tau must be > 0.5");
        }
        if (tau < 0.55 || tau > 1.2) {
            std::printf("LBM warning: tau out of stability range, simulation may be unstable.\n");
        }
    }; 

    /// Nx : Number of cells along the X axis
    /// Ny : Number of cells along the Y axis
    const int Nx, Ny;

    /// Re : Reynold number (viscosity related)
    /// u_lid : Velocity of the lid cavity
    const double Re, u_lid;

    const double w0 = 4.0/9.0;  // weight in (dx,dy)=(0,0)
    const double ws = 1.0/9.0;  // weight for adjacent points
    const double wd = 1.0/36.0; // diagonal weight

    // Arrays of the lattice weights and direction components
    const double wi[9] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};

    // direction numbering scheme:
    // 7 4 8
    // 3 0 1
    // 6 2 5
    const int dirx[9] = {0,1,0,-1, 0,1,-1,-1, 1};
    const int diry[9] = {0,0,1, 0,-1,1, 1,-1,-1};
    
    // The kinematic viscosity and the corresponding relaxation parameter
    const double nu;
    const double tau;
    //relaxaztion parameters
    const double taup_inv;
    const double taum_inv;

    // number of directions (considering the center point)
    const int ndir = 9; 
    // The fluid density
    const double rho0 = 1.0;

    // useful constants
    const double tauinv = 2.0/(6.0*nu+1.0);
    const double omtauinv = 1.0-tauinv;

    /**
    * Compute the equilibrium of the system.
    */
    void init_equilibrium();

    /**
    * Initialize the velocities and the density 
    * of the particle populations.
    */
    void init_lid_driven_cavity();

    /**
    * Apply the boundary conditions for
    * the lid cavity problem.
    */
    void apply_boundary_conditions();

    // TODO: ADD COMMENTS
    void update_stream_collide();

    // TODO: ADD COMMENTS
    void write_norms(std::ofstream& output_file);
    
    // TODO: ADD COMMENTS
    void write_bench_values(std::ofstream& output_file);

private:

    // Index position of a cell for a scalar defined vector
    inline size_t scalar_index(unsigned int x, unsigned int y) const {
        return Nx * y + x;
    }

    // Index position of a cell for a direction defined vector
    // This function is equal to: (Nx*Ny*d) + (Nx*y)+x
    // making d work as an offset.
    inline size_t field_index(unsigned int x, unsigned int y, unsigned int d) const {
        return Nx * (Ny * d + y) + x;
    }

    /**
     * Particle Populations
     */
    std::vector<double> f;
    /**
     * Temporary Particle Populations
     *
     * Used to store intermediate solutions
     */
    std::vector<double> f_tmp;
    /**
     * Velocitiy vector on the x axis
     */
    std::vector<double> ux;
    /**
     * Velocitiy vector on the x axis
     */
    std::vector<double> uy;
    /**
     * Vector of densities
     *
     * Used to store densities for benchmarking
     */
    std::vector<double> rho;
};

#endif