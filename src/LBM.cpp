#include "LBM.hpp"
#include <cmath>

//to do: implementation of functions for lattice boltzmann method

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

void LBM::stream(double *f_src, double* f_dst)
{
    for(unsigned int y = 0; y< Ny;++y)
    {
        for(unsigned int x = 0; x < Nx; ++x)
        {
            for(unsigned int i = 0; i < ndir; ++i)
            {
                // enforce periodicity: add Nx to ensure that value is positive
                // trova i vicini del vettore velocità u(x,y)
                unsigned int xmd = (Nx+x-dirx[i]) % Nx;
                unsigned int ymd = (Ny+y-diry[i]) % Ny;
                //per ogni peso del punto x,y assegno il peso associato del vicino xmd,ymd
                f_dst[field_index(x,y,i)] = f_src[field_index(xmd,ymd,i)];
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
                f[field_index(x,y,i)] =omtauinv*f[field_index(x,y,i)]+tauinv*feq;
            }
        }
    }
}

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



//to do: include necessary things

//initialize variables of the domain with class LBM



//write functions for: collision, macros quantity calculations, boundary conditions (bottom, top, right left)
void LBM::apply_boundary_conditions(double *u, double *v, double *r, double *f_src, double *f_dst) {

    // Left e Right boundaries
    for (int j = 1; j < Ny - 1; j++) {
        // Left boundary (x = 0)
        u[scalar_index(0, j)] = 0.0;
        v[scalar_index(0, j)] = 0.0;
        u[scalar_index(Nx-1, j)] = 0.0;
        v[scalar_index(Nx-1, j)] = 0.0;

        r[scalar_index(0,j)] = r[scalar_index(1,j)];
        r[scalar_index(Nx-1,j)] = r[scalar_index(Nx-2,j)];

        //DA RIVEDERE!!!
        for (int k = 0; k < ndir; k++) {
            f_dst[field_index(0, j, k)] = f_src[field_index(0, j, k)] + f_dst[field_index(1, j, k)] - f_src[field_index(1, j, k)];
            f_dst[field_index(Nx - 1, j, k)] = f_src[field_index(Nx - 1, j, k)] + f_dst[field_index(Nx - 2, j, k)] - f_src[field_index(Nx - 2, j, k)];
        }
    }

    // Top e Bottom boundaries
    for (int i = 0; i < Nx; i++) {
        // Top boundary (y = NY-1)
        u[scalar_index(i, Ny - 1)] = u_lid;
        v[scalar_index(i, Ny - 1)] = 0.0;
        u[scalar_index(i, 0)] = 0.0;
        v[scalar_index(i, 0)] = 0.0;

        r[scalar_index(i,Ny-1)] = r[scalar_index(i,Ny - 2)];
        r[scalar_index(i,0)] = r[scalar_index(i,1)];

        //DA RIVEDERE!!!
        for (int k = 0; k < ndir; k++) {
            f_dst[field_index(i, Ny - 1, k)] = f_src[field_index(i, Ny - 1, k)] + f_dst[field_index(i, Ny - 2, k)] - f_src[field_index(i, Ny - 2, k)];
            f_dst[field_index(i, 0, k)] = f_src[field_index(i, 0, k)] + f_dst[field_index(i, 1, k)] - f_src[field_index(i, 1, k)];
        }
    }
}

/*
void initializeGrid(){}
void generateVelocityField(){}
void simulateStep(){} //new f_temp during collision
double equilibriumDistribution(){}
double dotProduct(){}
void setBoundaryVelocities(){}
*/
