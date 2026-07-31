#pragma once

#include "lbm-2-lbm/problems/problem_base.hpp"

#include <omp.h>
#include <cmath>
#include <stdexcept>

namespace lbm {

template<enum ExecutionBackend exec_backend>
class LidCavity3D : public ProblemBase3D<exec_backend> {
public:
    using Base = ProblemBase3D<exec_backend>;
    using ExecutionContextT = typename Base::ExecutionContextT;

    void init(Grid<3>& grid, const double init_vel) const override {
        throw std::runtime_error("3D Lid cavity not yet implemented");
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
    void mpi_apply_boundary_conditions_3d(
        Grid<3>& grid, 
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context
    ) const {
        (void)grid;
        (void)part_stream;
        (void)context;
        throw std::runtime_error("3D MPI Lid cavity not yet implemented");
    };

    void cuda_apply_boundary_conditions_3d(
        Grid<3>& grid,
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context
    ) const {
        (void)grid;
        (void)part_stream;
        (void)context;
        throw std::runtime_error("3D Lid cavity not yet implemented");
    };

};

}
