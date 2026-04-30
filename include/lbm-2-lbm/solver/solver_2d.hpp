#pragma once

#include "lbm-2-lbm/solver/solver_base.hpp"
#include "lbm-2-lbm/types/common.hpp"

#include <omp.h>
#include <fstream>
#include <cmath>

namespace lbm {

template<class CollisionOperatorStrategy>
class MPISolver2D : public SolverBase2D<CollisionOperatorStrategy, ExecutionBackend::MPI> {
    using Base = SolverBase2D<CollisionOperatorStrategy, ExecutionBackend::MPI>;
public:
    MPISolver2D(
        Grid<2>& _grid,
        const double _init_vel,
        const Params<2, CollisionOperatorStrategy::TYPE> _params
    ) :
        Base(_grid, _init_vel, _params)
    {};

    ~MPISolver2D() = default;



    void init_equilibrium(std::vector<double>& part_stream) override {
        #pragma omp parallel for schedule(static) collapse(2)
        for(unsigned int y = 0; y < this->grid.Ny; ++y)
        {
            for(unsigned int x = 0; x < this->grid.Nx; ++x)
            {
                double r = this->grid.rho[this->grid.scalar_index(x,y)];
                double ux_val = this->grid.ux[this->grid.scalar_index(x,y)];
                double uy_val = this->grid.uy[this->grid.scalar_index(x,y)];
    
                #pragma omp simd
                for(unsigned int i = 0; i < D2Q9::ndir; ++i)
                {
                    double cidotu = 
                        D2Q9::dirx[i]*ux_val + D2Q9::diry[i]*uy_val;
                    // wi sono i coefficienti di peso del modello D2Q9
                    // Imposta f nella direzione i uguale alla distribuzione di equilibrio corrispondente a densità this->grid.rho e velocità (ux,uy).
                    part_stream[this->grid.field_index(x,y,i)] = 
                        D2Q9::wi[i]*r*(1.0 + 3.0*cidotu+4.5*cidotu*cidotu
                        -1.5*(ux_val*ux_val+uy_val*uy_val));
                }
            }
        }
    };


    void stream(
    std::vector<double>& ffrom,
    std::vector<double>& fto
    ) override {
        #pragma omp parallel for schedule(static) collapse(2)
        for (unsigned int y = 0; y < this->grid.Ny; ++y) {
            for (unsigned int x = 0; x < this->grid.Nx; ++x) {

                // ---- skip solid destination nodes ----
                // apply_bc (bounce-back) handles the fluid nodes
                // adjacent to the solid; solid internals stay at 0
                // (zeroed by update_collide on the previous step).
                if (this->grid.solid[this->grid.scalar_index(x, y)]) continue;

                #pragma omp unroll full
                for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
                    unsigned int xs = x - D2Q9::dirx[i];
                    unsigned int ys = y - D2Q9::diry[i];

                    if (xs < this->grid.Nx && ys < this->grid.Ny) {
                        fto[this->grid.field_index(x, y, i)] =
                            ffrom[this->grid.field_index(xs, ys, i)];
                    } else {
                        fto[this->grid.field_index(x, y, i)] = 0.0;
                    }
                }
            }
        }
    };


    void update_collide(
    std::vector<double>& ffrom,
    std::vector<double>& fto,
    bool save
    ) override {
    #pragma omp parallel for schedule(static) collapse(2)
        for (unsigned int y = 0; y < this->grid.Ny; ++y) {
            for (unsigned int x = 0; x < this->grid.Nx; ++x) {

                const unsigned int s_idx = this->grid.scalar_index(x, y);

                // ---- skip solid nodes ----
                // Zero out fto on solid so that after swap(f1,f2)
                // no stale populations leak into the next stream step.
                if (this->grid.solid[s_idx]) {
                    for (unsigned int i = 0; i < D2Q9::ndir; ++i)
                        fto[this->grid.field_index(x, y, i)] = 0.0;
                    continue;
                }

                double r      = 0.0;
                double ux_val = 0.0;
                double uy_val = 0.0;

                #pragma omp simd reduction(+:r,ux_val,uy_val)
                for (unsigned int i = 0; i < D2Q9::ndir; ++i) {
                    const double f_dist = fto[this->grid.field_index(x, y, i)];
                    r      += f_dist;
                    ux_val += D2Q9::dirx[i] * f_dist;
                    uy_val += D2Q9::diry[i] * f_dist;
                }

                ux_val /= r;
                uy_val /= r;

                if (save) {
                    this->grid.rho[s_idx] = r;
                    this->grid.ux [s_idx] = ux_val;
                    this->grid.uy [s_idx] = uy_val;
                }

                this->collision_strategy.apply(
                    this->grid,
                    r,
                    fto,
                    types::IndexPoint<2>(x, y),
                    types::ValuePoint<2>(ux_val, uy_val)
                );
            }
        }
    };
};

}
