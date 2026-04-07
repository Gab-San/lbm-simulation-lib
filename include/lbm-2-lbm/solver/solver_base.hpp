#pragma once

#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/core/metadata.hpp"
#include "lbm-2-lbm/types/core.hpp"
#include "lbm-2-lbm/core/velocity_sets.hpp"

#include <vector>

namespace lbm {

template<int dim, 
class VelocitySet, 
class CollisionOperatorStrategy, 
enum ExecutionBackend exec_backend>
class SolverBase {
protected:
    Grid<dim>& grid;
    const CollisionOperatorStrategy collision_strategy;

public:
    static constexpr types::CollisionModel CS_TYPE = CollisionOperatorStrategy::TYPE;
    static constexpr ExecutionBackend BACKEND_TYPE = exec_backend;

    SolverBase(
        Grid<dim>& _grid,
        const double _init_vel,
        const Params<dim, CollisionOperatorStrategy::TYPE> _params
    ): 
        grid(_grid),
        collision_strategy(_init_vel, _params)
    {};

    virtual ~SolverBase() = default;
    virtual void init_equilibrium(std::vector<double>& part_stream) = 0;
    virtual void stream(std::vector<double>& ffrom,
                        std::vector<double>& fto) = 0;
    virtual void update_collide(std::vector<double>& ffrom,
                                std::vector<double>& fto, bool save) = 0;
};

template<class CollisionOperatorStrategy, enum ExecutionBackend exec_backend>
class SolverBase2D : public SolverBase<2, D2Q9, CollisionOperatorStrategy, exec_backend> {
    using Base = SolverBase<2, D2Q9, CollisionOperatorStrategy, exec_backend>;

public:
    SolverBase2D(
        Grid<2>& _grid,
        const double _init_vel,
        const Params<2, CollisionOperatorStrategy::TYPE> _params
    ): 
        Base(_grid, _init_vel, _params)
    {};

    virtual ~SolverBase2D() = default;
    virtual void init_equilibrium(std::vector<double>& part_stream) override = 0;
    virtual void stream(std::vector<double>& ffrom, 
                        std::vector<double>& fto) override = 0;
    virtual void update_collide(std::vector<double>& ffrom,
                                std::vector<double>& fto, bool save) override = 0;
};

template<class VelocitySet, 
class CollisionOperatorStrategy, 
enum ExecutionBackend exec_backend>
class SolverBase3D : 
public SolverBase<3, VelocitySet, CollisionOperatorStrategy, exec_backend> {
    using Base = SolverBase<3, VelocitySet, CollisionOperatorStrategy, exec_backend>;

public:
    SolverBase3D(
        Grid<3>& _grid,
        const double _init_vel,
        const Params<3, CollisionOperatorStrategy::TYPE> _params
    ): 
        Base(_grid, _init_vel, _params)
    {};

    virtual ~SolverBase3D() = default;
    virtual void init_equilibrium(std::vector<double>& part_stream) override = 0;
    virtual void stream(std::vector<double>& ffrom, 
                        std::vector<double>& fto) override = 0;
    virtual void update_collide(std::vector<double>& ffrom,
                                std::vector<double>& fto, bool save) override = 0;
};

}
