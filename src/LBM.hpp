#ifndef LBM_HPP
#define LBM_HPP

#include "defs.hpp"

#include <stdexcept> // for runtime_error
#include <cstddef>   // for size_t
#include <cstdio>    // for printf
#include <array>
#include <vector>

/**
 * LBM
 *
 * \brief This class implements the lattice Boltzmann algorithm steps
 * considering a D2Q9 velocity set and the lid driven cavity problem.
 *
 * It contains the general parameters needed to implement the algorithm
 * and it is initialized with the parameters of the problem.
 *
 * It is flexible w.r.t. the size of the problem but
 * it is viable for the D2Q9 velocity set only.
 *
 * The steps used are:
 * - setup lid driven cavity problem
 * - initialize equilibrium
 * - stream, calculate densities and velocities & collide
 *
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
    LBM(int num_cells_x_, int num_cells_y_, double rey_num_, double u_lid_, lbm_lbm::CollisionOperator op_):
	op(op_),
        nu(u_lid_ * num_cells_y_ / rey_num_), 
        Nx(num_cells_x_), 
        Ny(num_cells_y_), 
        u_lid(u_lid_),
        Re(rey_num_),
	ux(num_cells_x_*num_cells_y_, 0.0),
	uy(num_cells_x_*num_cells_y_, 0.0),
	rho(num_cells_x_*num_cells_y_, rho0)
    {
        if (tau <= 0.5) {
            throw std::runtime_error("LBM error: tau must be > 0.5");
        }
        if (tau < 0.55 || tau > 1.2) {
            std::printf("LBM warning: tau out of stability range, simulation may be unstable.\n");
        }
    }; 

    /**
     * Number of direction of the D2Q9 set.
     */
    static const unsigned short int ndir = 9;
	
    /// Weight in (dx,dy)=(0,0)
    const double w0 = 4.0/9.0;

    /// Weight for adjacent points
    const double ws = 1.0/9.0;

    /// Diagonal weight
    const double wd = 1.0/36.0;

    /**
     * Direction weights map.
     *
     * The direction numbering scheme is: \n
     * ------ + x \n
     * |7 4 8 \n
     * |3 0 1 \n
     * |6 2 5 \n
     * +\n
     * y
     */
    const std::array<double, ndir> wi = {w0, ws, ws, ws, ws, wd, wd, wd, wd};

    /**
     * Array of directions in the x direction following the numbering scheme.
     *
     * The direction numbering scheme is: \n
     * ------ + x \n
     * |7 4 8 \n
     * |3 0 1 \n
     * |6 2 5 \n
     * \n
     * y
     */
    const std::array<int, ndir> dirx = {0,1,0,-1,0,1,-1,-1,1};

    /**
     * Array of direction in the y direction following the numbering scheme.
     *
     * The direction numbering scheme is: \n
     * ------ + x \n
     * |7 4 8 \n
     * |3 0 1 \n
     * |6 2 5 \n
     * \n
     * y
     */
    const std::array<int, ndir> diry = {0,0,1,0,-1,1,1,-1,-1};

    /// Collision operator
    const lbm_lbm::CollisionOperator op;

    /// Number of cells along the X axis
    const int Nx;
    /// Number of cells along the Y axis
    const int Ny; 
    ///Reynold number (viscosity related)
    const double Re; 
    /// Velocity of the lid cavity
    const double u_lid; 
    /// Kinematic viscosity
    const double nu; 
    /// Relaxation parameter (used in BGK)
    const double tau = 0.5 + 3.0 * nu; 

    /// Fluid default density
    const double rho0 = 1.0;

    /// 
    const double tauinv = 2.0/(6.0*nu+1.0);

    ///
    const double omtauinv = 1.0-tauinv; 

    /** 
     * Tau positive relaxation parameter.
     *
     * (considering expanded collision formula)
     *
     * \important TRT not supported
     *
     * N.B: This parameter has not been checked 
     */
    const double tau_pos_inv = 2.0/(6.0*nu+1.0);

    /** 
     * Tau negative relaxation parameter.
     *
     * (considering expanded collision formula)
     *
     * \important TRT not supported
     *
     * N.B: This parameter is not correct
     */
    const double tau_min_inv = 2.0/(6.0*nu+1.0); 

    /**
     * Compute the equilibrium of the system.
     *
     * @param f particle population
     */
    void init_equilibrium(std::vector<double>& f);
    /**
     * Initialize the velocities and the density 
     * of the particle populations
     */
    void init_lid_driven_cavity();
    /**
     * Apply the boundary conditions for
     * the lid cavity problem.
     *
     * @param f particle population
     */
    void apply_boundary_conditions(std::vector<double>& f);

    void update_stream_collide(std::vector<double>& f_src, std::vector<double>& f_dst, bool save);
    
    /**
     * Write the norm of the macroscopic velocities.
     *
     * @param output_file file stream on which to write
     */
    void write_norms(std::ofstream& output_file);

    /**
     * Write the values of the velocities in the
     * y direction on the horizontal axis of symmetry
     * of the square cavity.
     *
     * @param output_file file stream on which to write
     */
    void write_bench_data(std::ofstream& output_file);

private:

    const std::array<int, ndir> opp = {0, 3, 4, 1, 2, 6, 5, 7, 8};

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
