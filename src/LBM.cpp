#include "LBM.hpp"

#include <omp.h>
#include <fstream>
#include <cmath>


//initialize velocities and pressure fields for the lid driven cavity flow
void LBM::init_lid_driven_cavity()
{
    // set lid velocity at the top boundary
    for(unsigned int x = 0; x < Nx; ++x)
    {
        size_t sidx = scalar_index(x,Ny-1);
        ux[sidx] = u_lid;
    }
}

//compute equilibrium
//evolution of the steps (for cycle) -> UPDATE collision, macros and apply boundary conditions
void LBM::init_equilibrium(std::vector<double> & f)
{
    #pragma omp parallel for schedule(static) collapse(2)
    for(unsigned int y = 0; y< Ny; ++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            double r = rho[scalar_index(x,y)];
            double ux_val = ux[scalar_index(x,y)];
            double uy_val = uy[scalar_index(x,y)];

	    #pragma omp simd
            for(unsigned int i = 0; i < ndir; ++i)
            {
                double cidotu = dirx[i]*ux_val + diry[i]*uy_val;
                // wi sono i coefficienti di peso del modello D2Q9
                // Imposta f nella direzione i uguale alla distribuzione di equilibrio corrispondente a densità rho e velocità (ux,uy).
                f[field_index(x,y,i)] = wi[i]*r*(1.0 + 3.0*cidotu+4.5*cidotu*cidotu-1.5*(ux_val*ux_val+uy_val*uy_val));
            }
        }
    }
}

// Single function that performs streaming, boundary conditions, macro computation, and collision
void LBM::update_stream_collide(std::vector<double> &  f_src, std::vector<double> &  f_dst)
{
    #pragma omp parallel for schedule(static) collapse(2)
    for (unsigned int y = 0; y < Ny; ++y) {
        for (unsigned int x = 0; x < Nx; ++x) {

            // first apply streaming
	    #pragma omp unroll full
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

            double r = 0.0;
            double ux_val =  0.0;
            double uy_val =  0.0;

            #pragma omp simd reduction(+:r,ux_val,uy_val)
            for (unsigned int i = 0; i < ndir; ++i) {

                // precompute f_dist
                double f_dist = f_dst[field_index(x,y,i)];

                // calculate macroscopic variables
                r += f_dist;
                ux_val  += dirx[i]*f_dist;
                uy_val  += diry[i]*f_dist;
            }

            ux_val /= r;
            uy_val /= r;

            // assign computed macroscopic values
            const unsigned int s_idx = scalar_index(x,y);
            rho[s_idx] = r;
            ux[s_idx] = ux_val;
            uy[s_idx] = uy_val;

            // Collisione con SIMD
            const double u_sq = ux_val*ux_val + uy_val*uy_val;
            const double c1 = -1.5 * u_sq;
            
	    if(op == lbm_lbm::TRT) {
		#pragma omp simd
		for (unsigned int i = 0; i < ndir; ++i)
		{
		    int io = opp[i];
		    if (i > io) continue;   // evita doppio aggiornamento

		    double fi  = f_dst[field_index(x,y,i)];
		    double fio = f_dst[field_index(x,y,io)];

		    // parti simmetrica / antisimmetrica
		    double f_plus  = 0.5 * (fi + fio);
		    double f_minus = 0.5 * (fi - fio);

		    // prodotti scalari
		    double cidotu  = dirx[i]*ux_val  + diry[i]*uy_val;
		    double cidotuo = dirx[io]*ux_val + diry[io]*uy_val;

		    // equilibrio
		    double feq_i = wi[i] * r *
			(1.0 + 3.0*cidotu + 4.5*cidotu*cidotu + c1);

		    double feq_io = wi[io] * r *
			(1.0 + 3.0*cidotuo + 4.5*cidotuo*cidotuo + c1);

		    double feq_plus  = 0.5 * (feq_i + feq_io);
		    double feq_minus = 0.5 * (feq_i - feq_io);

		    // rilassamento TRT
		    double f_plus_new  = f_plus  - tau_pos_inv * (f_plus  - feq_plus);
		    double f_minus_new = f_minus - tau_min_inv * (f_minus - feq_minus);

		    // ricostruzione
		    f_dst[field_index(x,y,i )] = f_plus_new + f_minus_new;
		    f_dst[field_index(x,y,io)] = f_plus_new - f_minus_new;
		}
	    }

	    if (op == lbm_lbm::BGK) {
		#pragma omp simd
		for (unsigned int i = 0; i < ndir; ++i) {
		    // calculate dot product beetwen the velocity u(x,y) and the direction vector to its neighbour
		    double cidotu = dirx[i]*ux_val + diry[i]*uy_val;
		    // calculate equilibrium
		    double feq = wi[i] * r * (1.0 + 3.0*cidotu + 4.5*cidotu*cidotu + c1);           
		    // relax to equilibrium
		    f_dst[field_index(x,y,i)] = omtauinv * f_dst[field_index(x,y,i)] + tauinv * feq;
		}
	    }
        }
    }

}

void LBM::apply_boundary_conditions(std::vector<double> & f)
{
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            #pragma omp parallel for simd
            for (int y = 0; y < Ny; y++) {
		// LEFT wall (x = 0)
                const int x0 = 0;
                f[field_index(x0,y,1)] = f[field_index(x0,y,3)];
                f[field_index(x0,y,5)] = f[field_index(x0,y,7)];
                f[field_index(x0,y,8)] = f[field_index(x0,y,6)];
		// RIGHT wall (x = Nx-1)
                const int xn = Nx - 1;
                f[field_index(xn,y,3)] = f[field_index(xn,y,1)];
                f[field_index(xn,y,6)] = f[field_index(xn,y,8)];
                f[field_index(xn,y,7)] = f[field_index(xn,y,5)];
            }
        }

        #pragma omp section
        {
            #pragma omp parallel for
            for (int x = 0; x < Nx; x++) {
		// BOTTOM wall (y = 0) 
                const int y0 = 0;
                f[field_index(x,y0,2)] = f[field_index(x,y0,4)];
                f[field_index(x,y0,5)] = f[field_index(x,y0,7)];
                f[field_index(x,y0,6)] = f[field_index(x,y0,8)];
		// TOP wall (y = Ny-1)
                const int yn = Ny - 1;

                double rho =
                    (
                        f[field_index(x,yn,0)] +
                        f[field_index(x,yn,1)] +
                        f[field_index(x,yn,3)] +
                        2.0 * (
                            f[field_index(x,yn,2)] +
                            f[field_index(x,yn,5)] +
                            f[field_index(x,yn,6)]
                        )
                    ) / (1.0 + u_lid);

                f[field_index(x,yn,4)] = f[field_index(x,yn,2)];
                f[field_index(x,yn,7)] = f[field_index(x,yn,5)] - 0.5 * rho * u_lid;
                f[field_index(x,yn,8)] = f[field_index(x,yn,6)] + 0.5 * rho * u_lid;
            }
        }
    }
}

void LBM::write_norms(std::ofstream& output_file) {
	// saving euclidian norms of the vectors at each step 
	for (int j = 0; j < Ny; ++j) {
		for (int i = 0; i < Nx; ++i) {
		    double vx = ux[Nx * j + i];
		    double vy = uy[Nx * j + i];
		    double v  = sqrt(vx * vx + vy * vy);
		    output_file << v << "\n";
		}
	}

}

void LBM::write_bench_data(std::ofstream& output_file) {
	// Results for v-Velocity along Horizontal Line through Geometric Center of Cavity
	int j_center = Ny / 2;
	for (int i = 0; i < Nx; ++i) {
	    double v_center = uy[Nx * j_center + i];
	    output_file << v_center << "\n";
	}
}
