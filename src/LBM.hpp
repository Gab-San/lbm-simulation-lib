
#ifndef LBM_HPP
#define LBM_HPP

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

#include <vector>
#include <string>

// qui abbiamo la classe LBM con parti private e parti pubbliche
class LBM {

// pubblica: ciò che il codice esterno può usare
public:

    /**
     *  Constructor
     *  @param grid_x_ number of cells along the X axis
     *  @param grid_y_ number of cells along the Y axis
     *  @param rey_num_ : Reynold number (viscosity related)
     *  @param u_lid : Initial velocity of the lid cavity
     */
    LBM(int grid_x_, int grid_y_, double rey_num_, double u_lid_):
	Nx(grid_x_), Ny(grid_y_), Re(rey_num_), u_lid(u_lid_){}; 

    /// Nx : Number of cells along the X axis
    /// Ny : Number of cells along the Y axis
    int Nx, Ny;

    /// Re : Reynold number (viscosity related)
    /// u_lid : Initial velocity of the lid cavity
    double Re, u_lid;
    double omega;

    const double w0 = 4.0/9.0;  // zero weight
    const double ws = 1.0/9.0;  // adjacent weight
    const double wd = 1.0/36.0; // diagonal weight

    // Arrays of the lattice weights and direction components
    const double wi[9] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};

    // direction numbering scheme:
    // 6 2 5
    // 3 0 1
    // 7 4 8
    const int dirx[9] = {0,1,0,-1, 0,1,-1,-1, 1};
    const int diry[9] = {0,0,1, 0,-1,1, 1,-1,-1};
    // The kinematic viscosity and the corresponding relaxation parameter
    const double nu = 1.0/6.0; //u_lid * (Ny - 1) / Re
    const double tau = 3.0*nu+0.5;

    const int ndir = 9; // number of directions (considerando anche il centro)

    // The maximum flow speed
    const double u_max = 0.04;
    // The fluid density
    const double rho0 = 1.0;


    //qui sotto ora ci sono i cosiddetti "metodi di uso alto livello":

    // One full time step: collision + streaming + boundaries + macros
    //void step();

    // Run many steps
    //void run(int nSteps, int outputEvery, const std::string &outputFolder);

    // Write macroscopic fields at a given time
    void writeFields(const std::string &filename) const;

    // Extract centerline profiles for validation
    //void writeCenterlineProfiles(const std::string &filename) const;
    
    //void initializeLattice(); ==> uguale al taylor green
    void init_equilibrium(double *f, double *r, double *u, double *v);
    void stream(double *f_src, double* f_dst);
    void compute_rho_u(double *f, double *r, double *u, double *v);
    void collide(double *f, double *r, double *u, double *v);
    void taylor_green(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v);
    void taylor_green(unsigned int t, double *r,  double *u, double *v);
    void init_lid_driven_cavity(double *u, double *v, double *r);
    void apply_boundary_conditions(double *f);
    //void apply_boundary_conditions(double *u, double *v, double *r, double *f_src, double *f_dst);


// privata: ciò che il codice esterno non può vedere
private:
    // Distribution functions: f[q][i] with flattened 2D index
    std::vector<double> f;      // size: Nx * Ny * 9
    std::vector<double> fTemp;  // buffer for streaming

    // Macroscopic fields
    std::vector<double> rho;    // size: Nx * Ny
    std::vector<double> ux;     // size: Nx * Ny
    std::vector<double> uy;     // size: Nx * Ny

    // Internal helpers
    inline size_t scalar_index(unsigned int x, unsigned int y) const {
        return Nx * y + x;
    }

    inline size_t field_index(unsigned int x, unsigned int y, unsigned int d) const {
        return Nx * (Ny * d + y) + x;
    }
};

#endif
