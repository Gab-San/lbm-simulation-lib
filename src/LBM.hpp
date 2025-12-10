
// Versione precedente
//all headers go here
//#ifndef LBM_HPP
//#define LBM_HPP
//
//class LBM{
//    private:
//        //D2Q9 model
//        //velocity configuration
//        //equilibrium weights (vector for each direction)
//        //velocity imposed on the upper lid (U?)
//        //main domain variables (rho, f, u, F)
//        double *rho, *f, *u, *F;
//        //f = population, rho = scalar density, u = velocity,  F ? (credo forze esterne)
//        
//        //functions to calc: density, direction, velocity, field, f_equilibrium
//    public:
//        //dimension in x, y directions
//        //? boh?
//        LBM(unsigned int nx, unsigned int ny, double u_lid, double Re);
//
//};
//#endif

#pragma once
#include <vector>
#include <string>

// qui abbiamo la classe LBM con parti private e parti pubbliche
class LBM {

// pubblica: ciò che il codice esterno può usare
public:
    LBM(int nx, int ny, double Re, double u_lid);   // costruttore

    int Nx, Ny;
    double Re, uLid;
    double omega;

    const double w0 = 4.0/9.0;  // zero weight
    const double ws = 1.0/9.0;  // adjacent weight
    const double wd = 1.0/36.0; // diagonal weight
    // Arrays of the lattice weights and direction components
    const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
    const int dirx[] = {0,1,0,-1, 0,1,-1,-1, 1};
    const int diry[] = {0,0,1, 0,-1,1, 1,-1,-1};
    // The kinematic viscosity and the corresponding relaxation parameter
    const double nu = 1.0/6.0;
    const double tau = 3.0*nu+0.5;


    //qui sotto ora ci sono i cosiddetti "metodi di uso alto livello":

    // One full time step: collision + streaming + boundaries + macros
    //void step();

    // Run many steps
    //void run(int nSteps, int outputEvery, const std::string &outputFolder);

    // Write macroscopic fields at a given time
    void writeFields(const std::string &filename) const;

    // Extract centerline profiles for validation
    //void writeCenterlineProfiles(const std::string &filename) const;

    // Getters for debug or postprocessing
    int getNx() const { return Nx; }
    int getNy() const { return Ny; }


// privata: ciò che il codice esterno non può vedere
private:
    // Distribution functions: f[q][i] with flattened 2D index
    std::vector<double> f;      // size: Nx * Ny * 9
    std::vector<double> fTemp;  // buffer for streaming

    // Macroscopic fields
    std::vector<double> rho;    // size: Nx * Ny
    std::vector<double> ux;     // size: Nx * Ny
    std::vector<double> uy;     // size: Nx * Ny

    // Lattice weights and velocities
    double w[9];
    int ex[9];
    int ey[9];

    // Internal helpers
    int idx(int x, int y) const { return y * Nx + x; }
    int qIndex(int x, int y, int q) const { return (y * Nx + x) * 9 + q; }

    //void initializeLattice(); ==> uguale al taylor green
    void taylor_green();
    void collide();
    void stream();          
    void compute_rho_u();
    void applyBoundaries();
    void computeMacros();
    double equilibrium(int q, int i) const;
};
