#pragma once

#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/core/metadata.hpp"

#include "lbm-2-lbm/types/core.hpp"


namespace lbm {

template<int dim, 
enum ExecutionBackend exec_backend = ExecutionBackend::MPI>
class ProblemBase {
    static_assert(dim == 2 || dim == 3);
public:
    static const enum ExecutionBackend BACKEND_TYPE = exec_backend;
    using ExecutionContextT = ExecutionContext<exec_backend>;
    
    virtual void init(Grid<dim>& grid, const double init_vel) const = 0;

    virtual void apply_boundary_conditions(
        Grid<dim>& grid,
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context = ExecutionContextT{}
    ) const = 0;
};

template<enum ExecutionBackend exec_backend = ExecutionBackend::MPI>
class ProblemBase2D : public ProblemBase<2, exec_backend> {
public:
    using Base = ProblemBase<2, exec_backend>;
    using ExecutionContextT = typename Base::ExecutionContextT;
    
    virtual ~ProblemBase2D() = default;

    virtual void init(Grid<2>& grid, const double init_vel) const override = 0;
    
    virtual void apply_boundary_conditions(
        Grid<2>& grid,
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context = ExecutionContextT{}
    ) const override = 0;
};

template<enum ExecutionBackend exec_backend = ExecutionBackend::MPI>
class ProblemBase3D : public ProblemBase<3, exec_backend> {
public:
    using Base = ProblemBase<3, exec_backend>;
    using ExecutionContextT = typename Base::ExecutionContextT;

    virtual ~ProblemBase3D() = default;

    virtual void init(Grid<3>& grid, const double init_vel) const override = 0;
    
    virtual void apply_boundary_conditions(
        Grid<3>& grid,
        const double init_vel,
        std::vector<double>& part_stream,
        const ExecutionContextT& context = ExecutionContextT{}
    ) const override = 0;
};

}
