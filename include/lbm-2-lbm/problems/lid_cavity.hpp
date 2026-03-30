#pragma once

#include "lbm-2-lbm/problems/problem_base.hpp"

#include <omp.h>
#include <cmath>
#include <stdexcept>

namespace lbm {

template<int dim>
class LidCavityProblem : public ProblemBase<dim> {
public:
    LidCavityProblem(
        Grid<dim>& _grid,
        const double _init_vel  
    ):
	ProblemBase<dim>(_grid, _init_vel)
    {};

    void init() override {
	if constexpr (dim == 2)
	    init_2d();
	else
	    init_3d();
    };
    void apply_boundary_conditions(std::vector<double>& part_stream) override {
	if constexpr (dim == 2)
	    apply_boundary_conditions_2d(part_stream);
	else
	    apply_boundary_conditions_3d();
    };

private:
    void init_2d() {
	#pragma omp simd
        // set lid velocity at the top boundary
	for(unsigned int x = 0; x < this->grid.Nx; ++x) {
	    this->grid.ux[this->grid.scalar_index(x,this->grid.Ny-1)] = this->init_vel;
	}
    };
    void apply_boundary_conditions_2d(std::vector<double>& part_stream) {
        #pragma omp parallel sections
	{
        #pragma omp section
        {
	    const int x0 = 0;
            const int xn = this->grid.Nx - 1;
            #pragma omp parallel for simd
            for (int y = 0; y < this->grid.Ny; y++) {
		// LEFT wall (x = 0)
                part_stream[this->grid.field_index(x0,y,1)] = 
		    part_stream[this->grid.field_index(x0,y,3)];
                part_stream[this->grid.field_index(x0,y,5)] = 
		    part_stream[this->grid.field_index(x0,y,7)];
                part_stream[this->grid.field_index(x0,y,8)] = 
		    part_stream[this->grid.field_index(x0,y,6)];
		// RIGHT wall (x = this->grid.Nx-1)
                part_stream[this->grid.field_index(xn,y,3)] = 
		    part_stream[this->grid.field_index(xn,y,1)];
                part_stream[this->grid.field_index(xn,y,6)] = 
		    part_stream[this->grid.field_index(xn,y,8)];
                part_stream[this->grid.field_index(xn,y,7)] = 
		    part_stream[this->grid.field_index(xn,y,5)];
            }
        }

        #pragma omp section
        {
            const int y0 = 0;
            const int yn = this->grid.Ny - 1;
            #pragma omp parallel for
            for (int x = 0; x < this->grid.Nx; x++) {
		// BOTTOM wall (y = 0) 
                part_stream[this->grid.field_index(x,y0,2)] = 
			part_stream[this->grid.field_index(x,y0,4)];
                part_stream[this->grid.field_index(x,y0,5)] = 
			part_stream[this->grid.field_index(x,y0,7)];
                part_stream[this->grid.field_index(x,y0,6)] = 
			part_stream[this->grid.field_index(x,y0,8)];

		// TOP wall (y = this->grid.Ny-1)
                double rho =
                    (
                        part_stream[this->grid.field_index(x,yn,0)] +
                        part_stream[this->grid.field_index(x,yn,1)] +
                        part_stream[this->grid.field_index(x,yn,3)] +
                        2.0 * (
                            part_stream[this->grid.field_index(x,yn,2)] +
                            part_stream[this->grid.field_index(x,yn,5)] +
                            part_stream[this->grid.field_index(x,yn,6)]
                        )
                    ) / (1.0 + this->init_vel);

                part_stream[this->grid.field_index(x,yn,4)] = 
			part_stream[this->grid.field_index(x,yn,2)];
                part_stream[this->grid.field_index(x,yn,7)] = 
			part_stream[this->grid.field_index(x,yn,5)] - 
			0.5 * rho * this->init_vel;
                part_stream[this->grid.field_index(x,yn,8)] = 
			part_stream[this->grid.field_index(x,yn,6)] + 
			0.5 * rho * this->init_vel;
            }
        }
	}
    };

    void init_3d() {
	throw std::runtime_error("3D Lid cavity not yet implemented");
    };
    void apply_boundary_conditions_3d() {
	throw std::runtime_error("3D Lid cavity not yet implemented");
    };
};

}
