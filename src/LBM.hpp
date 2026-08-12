#ifndef LBM_HPP
#define LBM_HPP

#include "defs.h"

#include <stdexcept> // for runtime_error
#include <cstddef>   // for size_t
#include <cstdio>    // for printf
#include <array>
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
    LBM(int num_cells_x_, int num_cells_y_, double rey_num_, double u_lid_, lbm_lbm::CollisionOperator op_):
	op(op_),
        nu(u_lid_ * num_cells_y_ / rey_num_), 
        tau(0.5 + 3.0 * nu),
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

    static const unsigned short int ndir = 9;

    const double w0 = 4.0/9.0;  // weight in (dx,dy)=(0,0)
    const double ws = 1.0/9.0;  // weight for adjacent points
    const double wd = 1.0/36.0; // diagonal weight

    const std::array<double, ndir> wi = {w0, ws, ws, ws, ws, wd, wd, wd, wd};

    const std::array<int, ndir> dirx = {0,1,0,-1,0,1,-1,-1,1};
    const std::array<int, ndir> diry = {0,0,1,0,-1,1,1,-1,-1};

    const lbm_lbm::CollisionOperator op;

    const int Nx, Ny;
    const double Re, u_lid;
    const double nu;
    const double tau;

    const double rho0 = 1.0;

    const double tauinv = 2.0/(6.0*nu+1.0);
    const double omtauinv = 1.0-tauinv;

    // Separate relaxation times if needed
    const double tau_pos_inv = 2.0/(6.0*nu+1.0);
    const double tau_min_inv = 2.0/(6.0*nu+1.0);

    // Funzioni principali
    void init_equilibrium(std::vector<double>& f);
    void init_lid_driven_cavity();
    void apply_boundary_conditions(std::vector<double>& f);
    void update_stream_collide(std::vector<double>& f_src, std::vector<double>& f_dst, bool save);
    void write_norms(std::ofstream& output_file);
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
