#ifndef LBM_HPP
#define LBM_HPP

#include <stdexcept> // for runtime_error
#include <cstddef>   // for size_t

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
        nu(u_lid_ * (num_cells_y_ - 1) / rey_num_), 
        tau(0.5 + 3.0 * nu),
        Nx(num_cells_x_), 
        Ny(num_cells_y_), 
        u_lid(u_lid_),
        Re(rey_num_)
    {
        if (tau <= 0.5) {
            throw std::runtime_error("LBM error: tau must be > 0.5");
        }
        if (tau < 0.55 || tau > 1.2) {
            std::printf("LBM warning: tau out of stability range, simulation may be unstable.\n");
        }

        // TODO: Trovare la stabilità per il Mach number
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
    const double nu = u_lid * (Ny - 1) / Re;
    const double tau = 3.0 * nu + 0.5;
    


    const int ndir = 9; // number of directions (considering the center point)

    // The maximum flow speed
    const double u_max = 0.04;
    // The fluid density
    const double rho0 = 1.0;

    // FIXME: Are these methods needed?

    // Write macroscopic fields at a given time
    // void writeFields(const std::string &filename) const;

    // Extract centerline profiles for validation
    // void writeCenterlineProfiles(const std::string &filename) const;

    /**
    * Compute the equilibrium of the system.
    *
    * @param f particle population
    * @param r density vector
    * @param u x-axis velocity vector
    * @param v y-axis velocity vector
    */
    void init_equilibrium(double *f, double *r, double *u, double *v);

    /**
    * Compute the particle populations movement from
    * one cell to another distributing the velocities.
    *
    * It basically copies the distribution function in 
    * f_src to f_dst.
    *
    * @param f_src source particle population
    * @param f_dst destination particle population
    */
    void stream(double *f_src, double* f_dst);

    /**
    * Compute the density and the velocities for each cell.
    *
    * @param f particle population
    * @param r density vector
    * @param u x-axis velocity vector
    * @param v y-axis velocity vector
     */
    void compute_rho_u(double *f, double *r, double *u, double *v);

    /**
    * Evaluate the collision operator (in this case BGK)
    * on the particle populations using the density and velocity values
    * computed at the current step.
    *
    * @param f particle population
    * @param r density vector
    * @param u x-axis velocity vector
    * @param v y-axis velocity vector
    */
    void collide(double *f, double *r, double *u, double *v);

    // FIXME: Are these methods needed?
    void taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v);
    void taylor_green(unsigned int t, double *r,  double *u, double *v);
    //void apply_boundary_conditions(double *u, double *v, double *r, double *f_src, double *f_dst);

    /**
    * Initialize the velocities and the density 
    * of the particle populations
    *
    * @param r density vector
    * @param u x-axis velocity vector
    * @param v y-axis velocity vector
    */
    void init_lid_driven_cavity(double *u, double *v, double *r);

    /**
    * Apply the boundary conditions for
    * the lid cavity problem.
    *
    * @param f particle population
    */
    void apply_boundary_conditions(double *f);

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
};

#endif
