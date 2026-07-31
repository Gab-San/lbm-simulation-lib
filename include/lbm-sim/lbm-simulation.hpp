#ifndef __LBM_SIM_LBM_SIMULATION_HPP
#define __LBM_SIM_LBM_SIMULATION_HPP

// LBM SIM LIB
#include "lbm-sim/core/grid.hpp"

#include "lbm-sim/backend/metadata.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/solver-base.hpp"

// COLLISION DETECTION LIB
#include "collision-detection/core/types.hpp"

// C++ STANDARD LIB
#include <filesystem>
#include <fstream>
#include <iostream>

#include <cstdint>
#include <cstring>
#include <vector>

namespace lbm {

template <int dim, typename VelocitySet,
          enum CollisionModel cm_t = CollisionModel::BGK>
class LBMSimulation : public DataObservable {
private:
  static_assert(
      dim == VelocitySet::dim,
      "LBMSimulation: template parameter 'dim' must match VelocitySet::dim");

  Grid<dim> grid;
  const Params<dim, cm_t> params;

public:
  LBMSimulation(const CollisionDetection::types::DimPoint<dim> grid_dim_,
                const Params<dim, cm_t> params_)
      : grid(grid_dim_), params(params_) {};

  template <enum ExecutionBackend backend_t>
  void solve(SolverBase<dim, VelocitySet, cm_t, backend_t> &solver,
             const LidCavity2D &problem) {
    std::cout << "Initializing Simulation." << std::endl;

    std::vector<double> f1(grid.getArea() * VelocitySet::ndir, 0.0);
    std::vector<double> f2(grid.getArea() * VelocitySet::ndir, 0.0);

    // FIXME: cannot be initialized like this
    // define generic initialization

    // CollisionDetection::Segment<2> seg(
    //     CollisionDetection::types::Coordinate<2>(0, grid.size.y - 1),
    //     CollisionDetection::types::Coordinate<2>(grid.size.x - 1,
    //                                              grid.size.y - 1));
    //
    // problem.init(grid, params.init_vel, seg.getPerimeter());

    // FIXME: check that initialization + init_equilibrium suffices
    std::cout << "Problem Initialized." << std::endl;
    solver.init_equilibrium(grid, f1);
    std::cout << "Equilibrium Initialized." << std::endl;

    write_header(grid);

    solver.solve(grid, params, f1, f2);

    std::cout << "Finished Simulation." << std::endl;
  };

  void output(const char *filepath) {
    using namespace std::filesystem;

    path outpath(filepath);
    path parent = outpath.parent_path();

    if (!exists(parent))
      create_directories(parent);

    std::cout << "Opening " << outpath << std::endl;

    std::ofstream fout(outpath, std::ios::binary);

    if (!fout.is_open())
      std::cerr << "Failed to create file: " << outpath << std::endl;

    std::cout << "Writing..." << std::endl;

    std::vector<double> v_center(grid.size.x);
    int j_center = grid.size.y / 2;
    for (int i = 0; i < grid.size.x; ++i) {
      v_center[i] = grid.u[grid.size.x * j_center + i].dy;
    }

    fout.write(reinterpret_cast<const char *>(v_center.data()),
               v_center.size() * sizeof(double));

    std::cout << "Finished writing to " << outpath << std::endl;

    fout.close();
  };

private:
  void write_header(const Grid<dim> &grid) {
    std::vector<char> buf(sizeof(int32_t) * dim);

    if constexpr (dim == 2) {
      const int32_t nx32 = static_cast<int32_t>(grid.size.x);
      const int32_t ny32 = static_cast<int32_t>(grid.size.y);

      std::memcpy(buf.data(), &nx32, sizeof(int32_t));
      std::memcpy(buf.data() + sizeof(int32_t), &ny32, sizeof(int32_t));
    } else {
      const int32_t nx32 = static_cast<int32_t>(grid.size.x);
      const int32_t ny32 = static_cast<int32_t>(grid.size.y);
      const int32_t nz32 = static_cast<int32_t>(grid.size.z);

      std::memcpy(buf.data(), &nx32, sizeof(int32_t));
      std::memcpy(buf.data() + sizeof(int32_t), &ny32, sizeof(int32_t));
      std::memcpy(buf.data() + 2 * sizeof(int32_t), &nz32, sizeof(int32_t));
    }

    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm

#endif // __LBM_SIM_LBM_SIMULATION_HPP
