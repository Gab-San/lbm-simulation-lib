//step da fare:

//load file dat
/* std::vector<std::pair<double,double>> airfoilPoints;

//read airfoil.dat
void loadAirfoil(const std::string& filename) {
    std::ifstream file(filename);
    double x, y;

    while (file >> x >> y) {
        airfoilPoints.emplace_back(x, y);
    }
}

//riempi gli spazi
bool pointInPolygon(int x, int y,
    const std::vector<std::pair<int,int>>& poly) {

    bool inside = false;
    int n = poly.size();

    for (int i = 0, j = n - 1; i < n; j = i++) {
        auto [xi, yi] = poly[i];
        auto [xj, yj] = poly[j];

        bool intersect =
            ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (double)(yj - yi) + xi);

        if (intersect)
            inside = !inside;
    }

    return inside;
}

//segna le celle solide in PROBLEM BASE
Grid grid;

void markAirfoil() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            if (pointInPolygon(i, j, airfoilGridPoints)) {
                grid(i,j).setSolid(true);
            }
        }
    }
}

//in applybound cond
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {

        if (grid(i,j).isSolid()) {
            applyBounceBack(i, j);
        }
    }
}

//in probelm base sarà qualcosa tipo:
class AirfoilProblem : public ProblemBase {

public:

    void initialize() override {
        loadAirfoil("naca0012.dat");
        buildAirfoil();
        markAirfoil();
    }

    void applyBoundaryConditions() override {
        applyWalls();
        applyInlet();
        applyAirfoil();
    }

};