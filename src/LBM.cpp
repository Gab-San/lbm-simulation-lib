#include "LBM.hpp"
#include <cmath>


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

void LBM::apply_boundary_conditions(double* f)
{
    // LEFT wall (x = 0)
    for (int y = 0; y < Ny; y++) {
        int x = 0;
        f[field_index(x,y,1)] = f[field_index(x,y,3)];
        f[field_index(x,y,5)] = f[field_index(x,y,7)];
        f[field_index(x,y,8)] = f[field_index(x,y,6)];
    }

    // RIGHT wall (x = Nx-1)
    for (int y = 0; y < Ny; y++) {
        int x = Nx - 1;
        f[field_index(x,y,3)] = f[field_index(x,y,1)];
        f[field_index(x,y,6)] = f[field_index(x,y,8)];
        f[field_index(x,y,7)] = f[field_index(x,y,5)];
    }

    // BOTTOM wall (y = 0) 
    for (int x = 0; x < Nx; x++) {
        int y = 0;
        f[field_index(x,y,2)] = f[field_index(x,y,4)];
        f[field_index(x,y,5)] = f[field_index(x,y,7)];
        f[field_index(x,y,6)] = f[field_index(x,y,8)];
    }

    // TOP wall (y = Ny-1)
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

// FIXME: Are these functions needed?

void LBM::taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v)
{
    double kx = 2.0*M_PI/Nx;
    double ky = 2.0*M_PI/Ny;
    double td = 1.0/(nu*(kx*kx+ky*ky));
    double X = x+0.5;
    double Y = y+0.5;
    double ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    double uy = u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
    

    double P = -0.25*rho0*u_max*u_max*( (ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y) )*exp(-2.0*t/td);
    double rho = rho0+3.0*P;
    *r = rho;
    *u = ux;
    *v = uy;
}


void LBM::taylor_green(unsigned int t, double *r,  double *u, double *v)
{
    for(unsigned int y = 0; y < Ny; ++y)
        for(unsigned int x = 0; x < Nx; ++x)
        {
            size_t sidx = scalar_index(x,y);
            LBM::taylor_green(t,x,y,&r[sidx],&u[sidx],&v[sidx]);
        }
}
