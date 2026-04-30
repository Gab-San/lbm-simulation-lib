#pragma once

#include "lbm-2-lbm/types/common.hpp"

#include <cstddef>
#include <vector>

namespace lbm {

template<int dim>
struct Grid;

template<>
struct Grid<2> {
    /// Number of cells along the X axis
    const std::size_t Nx;
    /// Number of cells along the Y axis
    const std::size_t Ny;

    /// Velocity vector on the x axis
    std::vector<double> ux;
    /// Velocity vector on the y axis
    std::vector<double> uy;

    /// Vector of densities
    std::vector<double> rho;
    /// Fluid default density
    const double rho0 = 1.0;

    /// Solid obstacle mask — true means the node is inside a solid body.
    /// Initialised to false (all fluid). Problems that include solid
    /// obstacles (e.g. airfoils, cylinders) set this during init().
    std::vector<bool> solid;

    Grid(types::DimPoint<2> _grid_dim):
        Nx(_grid_dim.x), Ny(_grid_dim.y),
        ux (_grid_dim.x * _grid_dim.y, 0.0),
        uy (_grid_dim.x * _grid_dim.y, 0.0),
        rho(_grid_dim.x * _grid_dim.y, rho0),
        solid(_grid_dim.x * _grid_dim.y, false)
    {}

    /// Linear index for macroscopic fields (rho, ux, uy, solid)
    inline std::size_t scalar_index(std::size_t x, std::size_t y) const {
        return Nx * y + x;
    }

    /// Linear index into the distribution-function array f[].
    /// Layout is SoA: f is stored direction-major, i.e.
    ///   f[ dir * Nx*Ny  +  y * Nx  +  x ]
    /// which equals  Nx * (Ny * dir + y) + x
    inline std::size_t field_index(std::size_t x, std::size_t y, std::size_t dir) const {
        return Nx * (Ny * dir + y) + x;
    }

    inline std::size_t getArea() const {
        return Nx * Ny;
    }
};

template<>
struct Grid<3> {
    /// Number of cells along the X axis
    const std::size_t Nx;
    /// Number of cells along the Y axis
    const std::size_t Ny;
    /// Number of cells along the Z axis
    const std::size_t Nz;

    /// Velocity vector on the x axis
    std::vector<double> ux;
    /// Velocity vector on the y axis
    std::vector<double> uy;
    /// Velocity vector on the z axis
    std::vector<double> uz;

    /// Vector of densities
    std::vector<double> rho;
    /// Fluid default density
    const double rho0 = 1.0;

    /// Solid obstacle mask — true means the node is inside a solid body.
    std::vector<bool> solid;

    Grid(types::DimPoint<3> _grid_dim) :
        Nx(_grid_dim.x), Ny(_grid_dim.y), Nz(_grid_dim.z),
        ux (_grid_dim.x * _grid_dim.y * _grid_dim.z, 0.0),
        uy (_grid_dim.x * _grid_dim.y * _grid_dim.z, 0.0),
        uz (_grid_dim.x * _grid_dim.y * _grid_dim.z, 0.0),
        rho(_grid_dim.x * _grid_dim.y * _grid_dim.z, rho0),
        solid(_grid_dim.x * _grid_dim.y * _grid_dim.z, false)
    {}

    // TODO: CHECK THIS IMPLEMENTATION
    inline std::size_t scalar_index(
        std::size_t x, std::size_t y, std::size_t z
    ) const {
        return Nx * (Ny * z + y) + x;
    }

    inline std::size_t field_index(std::size_t x, std::size_t y, std::size_t z, std::size_t dir) const {
        return Nx * (Ny * (Nz * dir + z) + y) + x;
    }

    inline std::size_t getArea() const {
        return Nx * Ny * Nz;
    }
};

} // namespace lbm