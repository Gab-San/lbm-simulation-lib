#pragma once

#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/core/metadata.hpp"

#include "lbm-2-lbm/types/common.hpp"
#include "lbm-2-lbm/types/core.hpp"

#include "lbm-2-lbm/problems/problem_base.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <filesystem>

namespace lbm {

template<
int dim,
class VelocitySet,
class Solver,
enum types::CollisionModel cs_type = types::BGK,
enum ExecutionBackend exec_backend = ExecutionBackend::MPI
>
class LBMSimulation {
private:
    static_assert( dim == VelocitySet::dim,
        "LBMSimulation: template parameter 'dim' must match VelocitySet::dim" );
    static_assert( cs_type == Solver::CS_TYPE,
        "LBMSimulation: collision operator type must match Solver::CS_TYPE" );
    static_assert( exec_backend == Solver::BACKEND_TYPE,
        "LBMSimulation: backend must match solver backend");


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

    void solve(const ProblemBase<dim, exec_backend>& problem) {
        std::cout << "Initializing Simulation." << std::endl;

        Solver solver(grid, init_vel, params);
        std::vector<double> f1(grid.getArea() * VelocitySet::ndir, 0.0);
        std::vector<double> f2(grid.getArea() * VelocitySet::ndir, 0.0);

        problem.init(grid, init_vel);
        std::cout << "Problem Initialized." << std::endl;
        solver.init_equilibrium(f1);
        std::cout << "Equilibrium Initialized." << std::endl;

        #if listeners
        fw.write(GridSize);
        #endif

        if constexpr ( exec_backend == ExecutionBackend::MPI ) {
            for ( unsigned int iter = 0; iter < niters; iter++ ) {
                bool save = iter % nskips == 0;
                solver.stream(f1, f2);
                problem.apply_boundary_conditions(grid, init_vel, f2);
                solver.update_collide(f1, f2, save);
                std::swap(f1, f2);

                #if listeners
                if ( save ) fw.write(WriteNorms);
                #endif
            }
        } 
        else if constexpr (exec_backend == ExecutionBackend::CUDA)
        {
            // LAUNCH CUDA KERNEL
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
        create_directories(outpath.parent_path());

        std::ofstream fout(outpath, std::ios::binary);
        if (!fout.is_open()) {
            std::cerr << "Failed to create file: " << outpath << std::endl;
            return;
        }

        std::cout << "Writing to " << outpath << "..." << std::endl;

        const std::size_t Nx = grid.Nx;
        const std::size_t Ny = grid.Ny;

        // --- header: grid dimensions (2 x uint64) ---
        const uint64_t nx64 = Nx, ny64 = Ny;
        fout.write(reinterpret_cast<const char*>(&nx64), sizeof(uint64_t));
        fout.write(reinterpret_cast<const char*>(&ny64), sizeof(uint64_t));

        // --- solid mask (Nx*Ny x uint8, 1=solid 0=fluid) ---
        // Written separately so the reader can mask out solid nodes.
        std::vector<uint8_t> solid_out(Nx * Ny);
        for (std::size_t i = 0; i < Nx * Ny; ++i)
            solid_out[i] = grid.solid[i] ? 1u : 0u;
        fout.write(reinterpret_cast<const char*>(solid_out.data()),
                solid_out.size() * sizeof(uint8_t));

        // --- macroscopic fields: rho, ux, uy (Nx*Ny x double each) ---
        // Solid nodes are written as NaN so any accidental read is obvious.
        std::vector<double> rho_out(Nx * Ny);
        std::vector<double> ux_out (Nx * Ny);
        std::vector<double> uy_out (Nx * Ny);

        constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

        for (std::size_t y = 0; y < Ny; ++y) {
            for (std::size_t x = 0; x < Nx; ++x) {
                const std::size_t idx = grid.scalar_index(x, y);
                if (grid.solid[idx]) {
                    rho_out[idx] = NaN;
                    ux_out [idx] = NaN;
                    uy_out [idx] = NaN;
                } else {
                    rho_out[idx] = grid.rho[idx];
                    ux_out [idx] = grid.ux [idx];
                    uy_out [idx] = grid.uy [idx];
                }
            }
        }

        fout.write(reinterpret_cast<const char*>(rho_out.data()),
                rho_out.size() * sizeof(double));
        fout.write(reinterpret_cast<const char*>(ux_out.data()),
                ux_out.size()  * sizeof(double));
        fout.write(reinterpret_cast<const char*>(uy_out.data()),
                uy_out.size()  * sizeof(double));

        std::cout << "Done. Written "
                << Nx << "x" << Ny << " grid ("
                << Nx*Ny << " nodes)." << std::endl;

        fout.close();
    }
};

}
