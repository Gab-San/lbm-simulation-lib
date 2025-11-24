//all headers go here
#ifndef LBM_HPP
#define LBM_HPP

class LBM{
    private:
        //D2Q9 model
        //velocity configuration
        //equilibrium weights (vector for each direction)
        //velocity imposed on the upper lid (U?)
        //main domain variables (rho, f, u, F)
        double *rho, *f, *u, *F;
        //f= population, rho = scalar density, u = velocity,  F ? (credo forze esterne)
        
        //functions to calc: density, direction, velocity, field, f_equilibrium
    public:
        //dimension in x, y directions
        //? boh?
        LBM(unsigned int nx, unsigned int ny, double u_lid, double Re);

};
#endif