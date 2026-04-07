#pragma once

#include "lbm-2-lbm/problems/problem_base.hpp"

#include <omp.h>
#include <cmath>
#include <stdexcept>

namespace lbm {

template<enum ExecutionBackend exec_backend>
class LidCavity2D : public ProblemBase2D<exec_backend> {
public:
    using Base = ProblemBase2D<exec_backend>;
    using ExecutionContextT = typename Base::ExecutionContextT;

    virtual ~LidCavity2D() = default;

    void init(Grid<2>& grid, const double init_vel) const override {
         // the problem shall be set up by the CPU
        #pragma omp simd
        // set lid velocity at the top boundary
        for(unsigned int x = 0; x < grid.Nx; ++x) {
           grid.ux[grid.scalar_index(x,grid.Ny-1)] = init_vel;
        }
    };

    void apply_boundary_conditions(
        Grid<2>& grid, 
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context = ExecutionContextT{}
    ) const override {
        if constexpr (exec_backend == ExecutionBackend::MPI)
        {
            mpi_apply_boundary_conditions_2d(grid, init_vel, part_stream, context);
        }
        else if constexpr (exec_backend == ExecutionBackend::CUDA)
        {
            cuda_apply_boundary_conditions_2d(grid, init_vel, part_stream, context);
        }
        
    };

private:
    void mpi_apply_boundary_conditions_2d(
        Grid<2>& grid, 
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context
    ) const  {
        (void)context;
        #pragma omp parallel sections
        {
        #pragma omp section
        {
            const int x0 = 0;
            const int xn = grid.Nx - 1;
            #pragma omp parallel for simd
            for (int y = 0; y < grid.Ny; y++) {
                // LEFT wall (x = 0)
                part_stream[grid.field_index(x0,y,1)] = 
                    part_stream[grid.field_index(x0,y,3)];
                part_stream[grid.field_index(x0,y,5)] = 
                    part_stream[grid.field_index(x0,y,7)];
                part_stream[grid.field_index(x0,y,8)] = 
                    part_stream[grid.field_index(x0,y,6)];
                // RIGHT wall (x = grid.Nx-1)
                part_stream[grid.field_index(xn,y,3)] = 
                    part_stream[grid.field_index(xn,y,1)];
                part_stream[grid.field_index(xn,y,6)] = 
                    part_stream[grid.field_index(xn,y,8)];
                part_stream[grid.field_index(xn,y,7)] = 
                    part_stream[grid.field_index(xn,y,5)];
            }
        }

        #pragma omp section
        {
            const int y0 = 0;
            const int yn = grid.Ny - 1;
            #pragma omp parallel for
            for (int x = 0; x < grid.Nx; x++) {
                // BOTTOM wall (y = 0) 
                part_stream[grid.field_index(x,y0,2)] = 
                        part_stream[grid.field_index(x,y0,4)];
                part_stream[grid.field_index(x,y0,5)] = 
                        part_stream[grid.field_index(x,y0,7)];
                part_stream[grid.field_index(x,y0,6)] = 
                        part_stream[grid.field_index(x,y0,8)];

                // TOP wall (y = grid.Ny-1)
                double rho =
                    (
                        part_stream[grid.field_index(x,yn,0)] +
                        part_stream[grid.field_index(x,yn,1)] +
                        part_stream[grid.field_index(x,yn,3)] +
                        2.0 * (
                            part_stream[grid.field_index(x,yn,2)] +
                            part_stream[grid.field_index(x,yn,5)] +
                            part_stream[grid.field_index(x,yn,6)]
                        )
                    ) / (1.0 + init_vel);

                part_stream[grid.field_index(x,yn,4)] = 
                        part_stream[grid.field_index(x,yn,2)];
                part_stream[grid.field_index(x,yn,7)] = 
                        part_stream[grid.field_index(x,yn,5)] - 
                        0.5 * rho * init_vel;
                part_stream[grid.field_index(x,yn,8)] = 
                        part_stream[grid.field_index(x,yn,6)] + 
                        0.5 * rho * init_vel;
            }
        }
        }
    };

    void cuda_apply_boundary_conditions_2d(
        Grid<2>& grid,
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context
    ) const {
        (void)grid;
        (void)part_stream;
        (void)context;
        throw std::runtime_error("2D CUDA Lid cavity not yet implemented");
    };

};

}
