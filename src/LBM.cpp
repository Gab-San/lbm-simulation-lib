#include "LBM.hpp"
#include <cmath>

// enable OpenMP for parallelization
#include <omp.h>

const int LBM::opp[9] = {0, 3, 4, 1, 2, 6, 5, 7, 8};

//initialize velocities and pressure fields for the lid driven cavity flow
void LBM::init_lid_driven_cavity(double *u, double *v, double *r)
{
    for(unsigned int y = 0; y < Ny; ++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            size_t sidx = scalar_index(x,y);
	        u[sidx] = 0.0;
	        v[sidx] = 0.0;
            r[sidx] = rho0;
        }
    }

    // set lid velocity at the top boundary
    for(unsigned int x = 0; x < Nx; ++x)
    {
        size_t sidx = scalar_index(x,Ny-1);
        u[sidx] = u_lid;
    }
}

//compute equilibrium
//evolution of the steps (for cycle) -> UPDATE collision, macros and apply boundary conditions
void LBM::init_equilibrium(double *f, double *r, double *u, double *v)
{
    for(unsigned int y = 0; y< Ny; ++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux = u[scalar_index(x,y)];
            double uy = v[scalar_index(x,y)];
            for(unsigned int i = 0; i < ndir; ++i)
            {
                double cidotu = dirx[i]*ux + diry[i]*uy;
                // wi sono i coefficienti di peso del modello D2Q9
                // Imposta f nella direzione i uguale alla distribuzione di equilibrio corrispondente a densità rho e velocità (ux,uy).
                f[field_index(x,y,i)] = wi[i]*rho*(1.0 + 3.0*cidotu+4.5*cidotu*cidotu-1.5*(ux*ux+uy*uy));
            }
        }
    }
}

// Single function that performs streaming, boundary conditions, macro computation, and collision
void LBM::update_stream_collide(double * f_src, double * f_dst, double * r, double * u, double * v)  
{
    #pragma omp parallel for schedule(static) collapse(2)
    for (unsigned int y = 0; y < Ny; ++y) {
        for (unsigned int x = 0; x < Nx; ++x) {

            // first apply streaming

            for (unsigned int i = 0; i < ndir; ++i) {

                unsigned int xs = x - dirx[i];
                unsigned int ys = y - diry[i];

                // stream only if the source node is internal
                if (xs < Nx && ys < Ny) {
                    f_dst[field_index(x,y,i)] = f_src[field_index(xs,ys,i)];
                } else {
                    // otherwise, set the value to 0
                    f_dst[field_index(x,y,i)] = 0.0;                   
                }
            }
        }
    }

    // apply boundary conditions
    apply_boundary_conditions(f_dst);        
    
    // finally compute macroscopic variables and collide        
    #pragma omp parallel for schedule(static) collapse(2)
    for (unsigned int y = 0; y < Ny; ++y) {
        for (unsigned int x = 0; x < Nx; ++x) {

            double rho = 0.0;
            double ux =  0.0;
            double uy =  0.0;

            #pragma omp simd reduction(+:rho,ux,uy)
            for (unsigned int i = 0; i < ndir; ++i) {

                // precompute f_dist
                double f_dist = f_dst[field_index(x,y,i)];

                // calculate macroscopic variables
                rho += f_dist;
                ux  += dirx[i]*f_dist;
                uy  += diry[i]*f_dist;
            }

            ux /= rho;
            uy /= rho;

            // assign computed macroscopic values
            const unsigned int s_idx = scalar_index(x,y);
            r[s_idx] = rho;
            u[s_idx] = ux;
            v[s_idx] = uy;

            // Collisione con SIMD
            const double u_sq = ux*ux + uy*uy;
            const double c1 = -1.5 * u_sq;
            //TRT Collision
            #pragma omp simd
            for (unsigned int i = 0; i < ndir; ++i)
            {
                int io = opp[i];
            

                double fi  = f_dst[field_index(x,y,i)];
                double fio = f_dst[field_index(x,y,io)];

                // symmetric part and antusymmetric
                double f_plus  = 0.5 * (fi + fio);
                double f_minus = 0.5 * (fi - fio);

                                double cidotu  = dirx[i]*ux  + diry[i]*uy;
                double cidotuo = dirx[io]*ux + diry[io]*uy;

          
                double feq_i = wi[i] * rho *
                    (1.0 + 3.0*cidotu + 4.5*cidotu*cidotu + c1);

                double feq_io = wi[io] * rho *
                    (1.0 + 3.0*cidotuo + 4.5*cidotuo*cidotuo + c1);

                double feq_plus  = 0.5 * (feq_i + feq_io);
                double feq_minus = 0.5 * (feq_i - feq_io);

                // two relaxation parameters for stability and viscosità
                double f_plus_new  = f_plus  - taup_inv * (f_plus  - feq_plus);
                double f_minus_new = f_minus - taum_inv * (f_minus - feq_minus);

                
                f_dst[field_index(x,y,i )] = f_plus_new + f_minus_new;
                f_dst[field_index(x,y,io)] = f_plus_new - f_minus_new;
            }

            /*for (unsigned int i = 0; i < ndir; ++i) {
                // calculate dot product beetwen the velocity u(x,y) and the direction vector to its neighbour
                double cidotu = dirx[i]*ux + diry[i]*uy;
                // calculate equilibrium
                double feq = wi[i] * rho * (1.0 + 3.0*cidotu + 4.5*cidotu*cidotu + c1);           
                // relax to equilibrium
                f_dst[field_index(x,y,i)] = omtauinv * f_dst[field_index(x,y,i)] + tauinv * feq;
            }*/
        }
    }

}

void LBM::apply_boundary_conditions(double* f)
{
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // LEFT wall (x = 0)
            #pragma omp parallel for
            for (int y = 0; y < Ny; y++) {
                int x = 0;
                f[field_index(x,y,1)] = f[field_index(x,y,3)];
                f[field_index(x,y,5)] = f[field_index(x,y,7)];
                f[field_index(x,y,8)] = f[field_index(x,y,6)];
            }
        }
        
        #pragma omp section
        {
            // RIGHT wall (x = Nx-1)
            #pragma omp parallel for
            for (int y = 0; y < Ny; y++) {
                int x = Nx - 1;
                f[field_index(x,y,3)] = f[field_index(x,y,1)];
                f[field_index(x,y,6)] = f[field_index(x,y,8)];
                f[field_index(x,y,7)] = f[field_index(x,y,5)];
            }
        }

        #pragma omp section
        {
            // BOTTOM wall (y = 0) 
            #pragma omp parallel for
            for (int x = 0; x < Nx; x++) {
                int y = 0;
                f[field_index(x,y,2)] = f[field_index(x,y,4)];
                f[field_index(x,y,5)] = f[field_index(x,y,7)];
                f[field_index(x,y,6)] = f[field_index(x,y,8)];
            }
        }   

        #pragma omp section
        {
            // TOP wall (y = Ny-1)
            #pragma omp parallel for
            for (int x = 0; x < Nx; x++) {
                int y = Ny - 1;
                
                double rho =
                    (
                        f[field_index(x,y,0)] +
                        f[field_index(x,y,1)] +
                        f[field_index(x,y,3)] +
                        2.0 * (
                            f[field_index(x,y,2)] +
                            f[field_index(x,y,5)] +
                            f[field_index(x,y,6)]
                        )
                    ) / (1.0 + u_lid);

                f[field_index(x,y,4)] = f[field_index(x,y,2)];
                f[field_index(x,y,7)] = f[field_index(x,y,5)] - 0.5 * rho * u_lid;
                f[field_index(x,y,8)] = f[field_index(x,y,6)] + 0.5 * rho * u_lid;
            }
        }
    }
}

// Streaming step is not periodic, so we need to be careful at the boundaries
void LBM::stream(double *f_src, double *f_dst)
{
    // streaming vero e proprio
    for (unsigned int y = 0; y < Ny; ++y) {
        for (unsigned int x = 0; x < Nx; ++x) {
            for (unsigned int i = 0; i < ndir; ++i) {

                unsigned int xs = x - dirx[i];
                unsigned int ys = y - diry[i];

                // stream SOLO se il nodo sorgente è interno
                if (xs >= 0 && xs < Nx && ys >= 0 && ys < Ny) {
                    f_dst[field_index(x,y,i)] = f_src[field_index(xs,ys,i)];
                } else {
                    // altrimenti, imposto il valore a 0
                    f_dst[field_index(x,y,i)] = 0.0;
                }
            }
        }
    }
}

void LBM::compute_rho_u(double *f, double *r, double *u, double *v)
{
    for(unsigned int y = 0; y< Ny;++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;

            for(unsigned int i = 0; i < ndir; ++i)
            {
                rho += f[field_index(x,y,i)];
                ux += dirx[i]*f[field_index(x,y,i)];
                uy += diry[i]*f[field_index(x,y,i)];
            }

            r[scalar_index(x,y)] = rho;
            u[scalar_index(x,y)] = ux/rho;
            v[scalar_index(x,y)] = uy/rho;
        }
    }
}


// implementation of the BGK collision operator
void LBM::collide(double *f, double *r, double *u, double *v)
{
    // useful constants
    const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;
    for(unsigned int y = 0; y< Ny;++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux = u[scalar_index(x,y)];
            double uy = v[scalar_index(x,y)];
            // 1- 1/tau
            for(unsigned int i = 0; i < ndir; ++i)
            {
                // calculate dot product beetwen the velocity u(x,y) and the direction vector to its neighbour
                double cidotu = dirx[i]*ux + diry[i]*uy;
                // calculate equilibrium
                double feq = wi[i]*rho*(1.0 + 3.0*cidotu+4.5*cidotu*cidotu-1.5*(ux*ux+uy*uy));
                // relax to equilibrium
                f[field_index(x,y,i)] = omtauinv*f[field_index(x,y,i)]+tauinv*feq;
            }
        }
    }
}
