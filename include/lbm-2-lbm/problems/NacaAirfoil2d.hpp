#pragma once

#include "lbm-2-lbm/problems/problem_base.hpp"
#include "lbm-2-lbm/core/grid.hpp"
#include "lbm-2-lbm/core/velocity_sets.hpp"
#include "lbm-2-lbm/types/core.hpp"

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace lbm {

// ============================================================
//  AirfoilGeometry  —  loads and normalises a .dat profile
// ============================================================

struct AirfoilGeometry {
    std::vector<double> x_pts;
    std::vector<double> y_pts;

    // ----------------------------------------------------------
    // Parse a NACA .dat file.  Handles:
    //   • Selig format   — text title header, then x y pairs
    //   • Lednicer format — numeric "n_upper. n_lower." header
    //   • Plain           — bare x y pairs, no header
    // ----------------------------------------------------------
    static AirfoilGeometry from_dat(const std::string& filepath) {
        std::ifstream fin(filepath);
        if (!fin.is_open())
            throw std::runtime_error(
                "AirfoilGeometry: cannot open file: " + filepath);

        AirfoilGeometry geom;
        std::string line;
        bool header_consumed = false;

        while (std::getline(fin, line)) {
            const auto first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) continue;
            const std::string trimmed = line.substr(first);

            if (std::isalpha(static_cast<unsigned char>(trimmed[0]))) {
                header_consumed = true;
                continue;
            }

            std::istringstream ss(trimmed);
            double x, y;
            if (!(ss >> x >> y)) continue;

            // Lednicer numeric header: x > 1.5 (e.g. "32.  32.")
            if (!header_consumed && x > 1.5) {
                header_consumed = true;
                continue;
            }

            geom.x_pts.push_back(x);
            geom.y_pts.push_back(y);
        }

        if (geom.x_pts.empty())
            throw std::runtime_error(
                "AirfoilGeometry: no valid coordinates in " + filepath);

        geom.normalize();
        return geom;
    }

    // chord -> [0,1], vertical centre -> 0
    void normalize() {
        const double xmin = *std::min_element(x_pts.begin(), x_pts.end());
        const double xmax = *std::max_element(x_pts.begin(), x_pts.end());
        const double chord = xmax - xmin;
        if (chord < 1e-12)
            throw std::runtime_error("AirfoilGeometry: degenerate chord");

        const double ymin = *std::min_element(y_pts.begin(), y_pts.end());
        const double ymax = *std::max_element(y_pts.begin(), y_pts.end());
        const double ymid = 0.5 * (ymin + ymax);

        for (std::size_t i = 0; i < x_pts.size(); ++i) {
            x_pts[i] = (x_pts[i] - xmin) / chord;
            y_pts[i] = (y_pts[i] - ymid) / chord;
        }
    }
};

// ============================================================
//  AirfoilRasterizer  —  marks grid.solid[] nodes
//
//  Parameters (as fractions of lattice dimensions):
//    chord_frac  : chord / Nx             (e.g. 0.30)
//    le_x_frac   : leading-edge x / Nx   (e.g. 0.25)
//    cy_frac     : airfoil centre y / Ny  (e.g. 0.50)
//    aoa_deg     : angle of attack [deg]  (e.g. 5.0)
// ============================================================

struct AirfoilRasterizer {
    static void rasterize(
        Grid<2>&               grid,
        const AirfoilGeometry& geom,
        double chord_frac = 0.30,
        double le_x_frac  = 0.25,
        double cy_frac    = 0.50,
        double aoa_deg    = 0.0
    ) {
        const int Nx = static_cast<int>(grid.Nx);
        const int Ny = static_cast<int>(grid.Ny);

        const double chord_px = chord_frac * Nx;
        const double le_x     = le_x_frac  * Nx;
        const double cy       = cy_frac     * Ny;

        // Rotation matrix for AoA (rotate profile, keep flow horizontal)
        const double aoa_rad = aoa_deg * M_PI / 180.0;
        const double cos_a   = std::cos(-aoa_rad);
        const double sin_a   = std::sin(-aoa_rad);

        // Rotate all profile points around mid-chord (x=0.5, y=0)
        const std::size_t N = geom.x_pts.size();
        std::vector<double> rx(N), ry(N);
        for (std::size_t i = 0; i < N; ++i) {
            const double px = geom.x_pts[i] - 0.5;
            const double py = geom.y_pts[i];
            rx[i] = cos_a * px - sin_a * py + 0.5;
            ry[i] = sin_a * px + cos_a * py;
        }

        // Map normalised coords -> lattice pixel coords
        auto to_px = [&](double nx) -> double { return le_x + nx * chord_px; };
        auto to_py = [&](double ny) -> double { return cy  + ny * chord_px; };

        // Per-column y-extent of the airfoil polygon (scan all edges)
        std::vector<double> y_lo(Nx,  1e30);
        std::vector<double> y_hi(Nx, -1e30);

        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t j = (i + 1) % N;

            const double x0 = to_px(rx[i]),  y0 = to_py(ry[i]);
            const double x1 = to_px(rx[j]),  y1 = to_py(ry[j]);

            const int xi0 = std::max(0,    static_cast<int>(std::floor(std::min(x0,x1))));
            const int xi1 = std::min(Nx-1, static_cast<int>(std::ceil (std::max(x0,x1))));

            for (int xi = xi0; xi <= xi1; ++xi) {
                const double dx = x1 - x0;
                const double t  = (std::abs(dx) < 1e-12)
                                  ? 0.5
                                  : std::clamp((xi - x0) / dx, 0.0, 1.0);
                const double yi = y0 + t * (y1 - y0);
                y_lo[xi] = std::min(y_lo[xi], yi);
                y_hi[xi] = std::max(y_hi[xi], yi);
            }
        }

        // Fill solid mask using grid.scalar_index(x, y)
        int solid_count = 0;
        for (int xi = 0; xi < Nx; ++xi) {
            if (y_lo[xi] > y_hi[xi]) continue;
            const int yi0 = std::max(0,    static_cast<int>(std::floor(y_lo[xi])));
            const int yi1 = std::min(Ny-1, static_cast<int>(std::ceil (y_hi[xi])));
            for (int yi = yi0; yi <= yi1; ++yi) {
                grid.solid[grid.scalar_index(xi, yi)] = true;
                ++solid_count;
            }
        }

        std::cout << "[AirfoilRasterizer] "
                  << solid_count << " solid nodes marked"
                  << "  (chord=" << static_cast<int>(chord_px) << "px"
                  << ", AoA=" << aoa_deg << "deg"
                  << ", grid=" << Nx << "x" << Ny << ")\n";
    }
};

// ============================================================
//  NACAirfoil2D  —  ProblemBase2D implementation
//
//  D2Q9 standard ordering:
//
//   dir │  0   1   2   3   4   5   6   7   8
//   ────┼────────────────────────────────────
//   cx  │  0  +1   0  -1   0  +1  -1  -1  +1
//   cy  │  0   0  +1   0  -1  +1  +1  -1  -1
//   opp │  0   3   4   1   2   7   8   5   6
//
//  f stored SoA (direction-major):
//     f[ grid.field_index(x, y, dir) ]
//       = Nx*(Ny*dir + y) + x
// ============================================================

template<enum ExecutionBackend exec_backend = ExecutionBackend::MPI>
class NACAirfoil2D : public ProblemBase2D<exec_backend> {
public:
    using Base              = ProblemBase2D<exec_backend>;
    using ExecutionContextT = typename Base::ExecutionContextT;

    static constexpr int Q = 9;
    static constexpr std::array<int, 9> CX  = { 0, 1, 0,-1, 0, 1,-1,-1, 1};
    static constexpr std::array<int, 9> CY  = { 0, 0, 1, 0,-1, 1, 1,-1,-1};
    static constexpr std::array<int, 9> OPP = { 0, 3, 4, 1, 2, 7, 8, 5, 6};

    // ---- configuration ----
    std::string dat_filepath;
    double chord_frac   = 0.30;
    double le_x_frac    = 0.25;
    double cy_frac      = 0.50;
    double aoa_deg      = 0.0;
    bool   open_channel = true; // true  -> free-slip top/bottom
                                // false -> periodic (stream wraps Ny)

    explicit NACAirfoil2D(
        const std::string& filepath,
        double chord   = 0.30,
        double le_x    = 0.25,
        double cy      = 0.50,
        double aoa     = 0.0,
        bool   open_ch = true
    ) :
        dat_filepath(filepath),
        chord_frac(chord), le_x_frac(le_x),
        cy_frac(cy), aoa_deg(aoa),
        open_channel(open_ch)
    {}

    virtual ~NACAirfoil2D() = default;

    // ----------------------------------------------------------
    //  init  —  load geometry, rasterize, initialise flow field
    // ----------------------------------------------------------
    void init(Grid<2>& grid, const double init_vel) const override {
        const int Nx = static_cast<int>(grid.Nx);
        const int Ny = static_cast<int>(grid.Ny);

        const auto geom = AirfoilGeometry::from_dat(dat_filepath);
        AirfoilRasterizer::rasterize(grid, geom,
                                     chord_frac, le_x_frac,
                                     cy_frac, aoa_deg);

        // Uniform free-stream: rho=1, u=(init_vel,0) for fluid nodes
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                const std::size_t idx = grid.scalar_index(x, y);
                grid.rho[idx] = 1.0;
                if (grid.solid[idx]) {
                    grid.ux[idx] = 0.0;
                    grid.uy[idx] = 0.0;
                } else {
                    grid.ux[idx] = init_vel;
                    grid.uy[idx] = 0.0;
                }
            }
        }
    }

    // ----------------------------------------------------------
    //  apply_boundary_conditions
    //
    //  Called every step, after stream(), before update_collide().
    //  f = post-streaming distribution (f2 in main loop).
    //  SoA access: f[ grid.field_index(x, y, dir) ]
    // ----------------------------------------------------------
    void apply_boundary_conditions(
        Grid<2>&              grid,
        const double          init_vel,
        std::vector<double>&  f,
        const ExecutionContextT& /*ctx*/ = ExecutionContextT{}
    ) const override {
        const int Nx = static_cast<int>(grid.Nx);
        const int Ny = static_cast<int>(grid.Ny);

        // SoA accessor  f(x, y, d)
        auto F = [&](int x, int y, int d) -> double& {
            return f[grid.field_index(x, y, d)];
        };

        // ---- 1. Airfoil: full-way bounce-back ----------------
        //
        // For each fluid node adjacent to a solid, populations
        // that tried to stream into the solid are sent back:
        //   F(x, y, opp[d]) <- F(x, y, d)
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                if (grid.solid[grid.scalar_index(x, y)]) continue;

                for (int d = 1; d < Q; ++d) {
                    const int nx = x + CX[d];
                    const int ny = y + CY[d];
                    if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny) continue;
                    if (!grid.solid[grid.scalar_index(nx, ny)]) continue;

                    F(x, y, OPP[d]) = F(x, y, d);
                }
            }
        }

        // ---- 2. Inlet: Zou-He velocity BC  (x = 0) ----------
        //
        // Imposed: ux = init_vel, uy = 0
        // Known dirs at left wall: 3, 6, 7
        // Unknown dirs:            1, 5, 8
        //
        // Zou & He (1997) left-wall relations:
        //   rho = (f0+f2+f4 + 2*(f3+f6+f7)) / (1 - ux)
        //   f1  = f3 + (2/3)*rho*ux
        //   f5  = f7 - (1/2)*(f2-f4) + (1/6)*rho*ux + (1/2)*rho*uy
        //   f8  = f6 + (1/2)*(f2-f4) + (1/6)*rho*ux - (1/2)*rho*uy
        for (int y = 1; y < Ny - 1; ++y) {
            const int x = 0;
            if (grid.solid[grid.scalar_index(x, y)]) continue;

            const double ux_w = init_vel;
            const double uy_w = 0.0;

            const double rho_w =
                (F(x,y,0) + F(x,y,2) + F(x,y,4)
               + 2.0 * (F(x,y,3) + F(x,y,6) + F(x,y,7)))
               / (1.0 - ux_w);

            F(x,y,1) = F(x,y,3) + (2.0/3.0)*rho_w*ux_w;
            F(x,y,5) = F(x,y,7)
                     - 0.5*(F(x,y,2) - F(x,y,4))
                     + (1.0/6.0)*rho_w*ux_w
                     + 0.5*rho_w*uy_w;
            F(x,y,8) = F(x,y,6)
                     + 0.5*(F(x,y,2) - F(x,y,4))
                     + (1.0/6.0)*rho_w*ux_w
                     - 0.5*rho_w*uy_w;

            grid.rho[grid.scalar_index(x,y)] = rho_w;
            grid.ux [grid.scalar_index(x,y)] = ux_w;
            grid.uy [grid.scalar_index(x,y)] = uy_w;
        }

        // ---- 3. Outlet: zero-gradient BC  (x = Nx-1) --------
        //
        // Copy all populations from last interior column (x = Nx-2).
        for (int y = 1; y < Ny - 1; ++y) {
            const int x = Nx - 1;
            if (grid.solid[grid.scalar_index(x, y)]) continue;

            double rho = 0.0, ux = 0.0, uy = 0.0;
            for (int d = 0; d < Q; ++d) {
                F(x, y, d) = F(x-1, y, d);
                rho += F(x, y, d);
                ux  += F(x, y, d) * CX[d];
                uy  += F(x, y, d) * CY[d];
            }
            grid.rho[grid.scalar_index(x,y)] = rho;
            grid.ux [grid.scalar_index(x,y)] = ux / rho;
            grid.uy [grid.scalar_index(x,y)] = uy / rho;
        }

        // ---- 4. Top / Bottom walls ---------------------------
        if (open_channel) {
            // Specular (free-slip) reflection:
            //
            // Bottom wall (y=0) — incoming dirs with cy<0: 4,7,8
            //   reflect to their cy-mirrored counterparts:  2,5,6
            for (int x = 0; x < Nx; ++x) {
                if (grid.solid[grid.scalar_index(x, 0)]) continue;
                F(x, 0, 2) = F(x, 0, 4);   // dir 4 -> dir 2
                F(x, 0, 5) = F(x, 0, 8);   // dir 8 -> dir 5
                F(x, 0, 6) = F(x, 0, 7);   // dir 7 -> dir 6
            }

            // Top wall (y=Ny-1) — incoming dirs with cy>0: 2,5,6
            //   reflect to: 4,8,7
            for (int x = 0; x < Nx; ++x) {
                if (grid.solid[grid.scalar_index(x, Ny-1)]) continue;
                F(x, Ny-1, 4) = F(x, Ny-1, 2);
                F(x, Ny-1, 8) = F(x, Ny-1, 5);
                F(x, Ny-1, 7) = F(x, Ny-1, 6);
            }
        }
        // else: periodic — stream() wraps y indices mod Ny

        // ---- 5. Domain corners: copy from nearest interior ---
        auto fix_corner = [&](int cx, int cy, int sx, int sy) {
            if (grid.solid[grid.scalar_index(cx, cy)]) return;
            for (int d = 0; d < Q; ++d)
                F(cx, cy, d) = F(sx, sy, d);
        };
        fix_corner(0,    0,    1,    1   );
        fix_corner(Nx-1, 0,    Nx-2, 1   );
        fix_corner(0,    Ny-1, 1,    Ny-2);
        fix_corner(Nx-1, Ny-1, Nx-2, Ny-2);
    }
};

} // namespace lbm