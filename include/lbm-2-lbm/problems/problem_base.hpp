#pragma once

#include "lbm-2-lbm/core/grid.hpp"

namespace lbm {

template<int dim>
class ProblemBase {
    static_assert(dim == 2 || dim == 3);
protected:
    Grid<dim>& grid;
    const double init_vel;

public:
    ProblemBase(
        Grid<dim>& _grid,
	const double _init_vel
    ):
        grid(_grid),
	init_vel(_init_vel)
    {};

    virtual void init() = 0;
    virtual void apply_boundary_conditions(std::vector<double>& part_stream) = 0;
};

}
