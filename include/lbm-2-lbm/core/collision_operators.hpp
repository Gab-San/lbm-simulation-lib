#pragma once

#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/core/metadata.hpp"

#include "lbm-2-lbm/types/common.hpp"
#include "lbm-2-lbm/types/core.hpp"

#include <vector>
#include <iostream>
#include <stdexcept>

namespace lbm {

template<int dim, class VelocitySet>
class CollisionStrategyBase {
    virtual void apply (
        Grid<dim>& grid, double r,
        std::vector<double>& fto, 
        types::IndexPoint<dim> p, 
        types::ValuePoint<dim> val
    ) const = 0;
};

template<int dim, class VelocitySet>
class BGKCollisionStrategy : public CollisionStrategyBase<dim, VelocitySet> {
    const double reyn_num;
    const double nu;
    const double tauinv;
    const double omtauinv;

public:
    static constexpr types::CollisionModel TYPE = types::BGK;

    BGKCollisionStrategy(const double _init_vel, const Params<dim, types::BGK>& params) :
        reyn_num(params.reyn_num),
        nu([&]() -> double {
            if constexpr (dim == 2) {
                return _init_vel * params.num_cells_y / params.reyn_num;
            } else {
                throw std::runtime_error("BGK: 3D not implemented yet!");
            }
        }()),
        tauinv(2.0 / (6.0 * nu + 1.0)),
        omtauinv(1.0 - tauinv)
    {
        double tau = 0.5 + 3.0 * nu;
        if (tau <= 0.5)
            throw std::runtime_error("LBM error: tau must be > 0.5");

        if (tau < 0.55 || tau > 1.2)
            std::cerr <<
            "LBM warning: tau out of stability range, simulation may be unstable." 
            << std::endl;
    };

    void apply(
		Grid<dim>& grid, double r,
		std::vector<double>& fto, 
		types::IndexPoint<dim> p, 
		types::ValuePoint<dim> val
    ) const override {
		if constexpr (dim == 2)
		    apply_2d(grid, r, fto, p, val);
		else
		    apply_3d(grid, r, fto, p, val);
    }

private:

    void apply_2d(
		Grid<2>& grid, double r,
		std::vector<double>& fto, 
		types::IndexPoint<2> p, 
		types::ValuePoint<2> val
    ) const {
		auto [x, y] = p;
		auto [ux_val, uy_val] = val;
		// Collisione con SIMD
		const double u_sq = ux_val*ux_val + uy_val*uy_val;
		// TODO: RENAME VARIABLE
		const double c1 = -1.5 * u_sq;

		#pragma omp simd
		for (unsigned int i = 0; i < VelocitySet::ndir; ++i) {
		    // calculate dot product beetwen the velocity u(x,y) 
		    // and the direction vector to its neighbour
		    double cidotu = VelocitySet::dirx[i]*ux_val + VelocitySet::diry[i]*uy_val;
		    // calculate equilibrium
		    double feq = VelocitySet::wi[i] * r * 
			(1.0 + 3.0*cidotu + 4.5*cidotu*cidotu + c1);
		    // relax to equilibrium
		    fto[grid.field_index(x,y,i)] = omtauinv * fto[grid.field_index(x,y,i)] + 
			tauinv * feq;
		}

    }

    void apply_3d(
		Grid<3>& grid, double r,
		std::vector<double>& fto, 
		types::IndexPoint<3> p, 
		types::ValuePoint<3> val
    ) const {
		auto [x,y,z] = p;
		throw std::runtime_error("3D version not implemented");
    }
};

template<int dim, class VelocitySet>
class TRTCollisionStrategy : public CollisionStrategyBase<dim, VelocitySet> {
public:
    static constexpr types::CollisionModel TYPE = types::CollisionModel::TRT;

    void apply(
		Grid<dim>& grid, double r,
		std::vector<double>& fto, 
		types::IndexPoint<dim> p, 
		types::ValuePoint<dim> val
    ) const override {
	/* TRT math */ 
    }

private:
    void apply_2d() const {};
    void apply_3d() const {};

};

}
