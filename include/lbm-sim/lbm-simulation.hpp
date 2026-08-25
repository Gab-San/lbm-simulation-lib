#ifndef __LBM_SIM_LBM_SIMULATION_HPP
#define __LBM_SIM_LBM_SIMULATION_HPP

#include "lbm-sim/core/grid.hpp"
#include "lbm-sim/core/types.hpp"

#include "lbm-sim/backend.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/problems/problem_2d.hpp"

#include "lbm-sim/solver/solver-base.hpp"

// C++ STANDARD LIB
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

namespace lbm {

template <unsigned short int dim, typename VelocitySet,
          enum CollisionModel cm_t = CollisionModel::BGK>
class LBMSimulation : public DataObservable {
private:
  static_assert(
      dim == VelocitySet::dim,
      "LBMSimulation: template parameter 'dim' must match VelocitySet::dim");

  Lattice<dim> lattice;
  const CollisionParams<dim, cm_t> params;

public:
  LBMSimulation(const types::DimPoint<dim> grid_dim_,
                types::boundary_mask_t boundary_mask_,
                const CollisionParams<dim, cm_t> params_, const double pin = 0,
                const double pout = 0)
      : lattice(grid_dim_, std::move(boundary_mask_), pin, pout),
        params(params_) {};

  template <enum ExecutionBackend backend_t>
  void solve(SolverBase<dim, VelocitySet, cm_t, backend_t> &solver,
             const LidCavity2D &problem) {

    // FIXME: cannot be initialized like this
    // define generic initialization

    // Segment<2> seg(
    //     types::Coordinate<2>(0, grid.size.y - 1),
    //     types::Coordinate<2>(grid.size.x - 1,
    //                                              grid.size.y - 1));
    //
    // problem.init(grid, params.init_vel, seg.getPerimeter());

    // FIXME: check that initialization + init_equilibrium suffices

    std::cout << "Initialized Simulation." << std::endl;
    write_header(lattice.grid);

    solver.solve(lattice, params);

    std::cout << "Finished Simulation." << std::endl;
  };

  void output(const char *filepath,
              std::function<std::vector<double>(const Lattice<dim> &)>
                  extract_profile) {
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

    std::string header = "%%profile " + collision_model_to_string(cm_t) + " " +
                         std::to_string(lattice.grid.size.y) + "\n";

    fout.write(header.data(), header.size());

    std::vector<double> profile = extract_profile(lattice);

    fout.write(reinterpret_cast<const char *>(profile.data()),
               profile.size() * sizeof(double));

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
