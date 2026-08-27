#ifndef __LBM_SIM_LBM_SIMULATION_HPP
#define __LBM_SIM_LBM_SIMULATION_HPP

#include "lbm-sim/types/common.hpp"

#include "lbm-sim/core/grid.hpp"

#include "lbm-sim/backend.hpp"
#include "lbm-sim/collision-operators/metadata.hpp"

#include "lbm-sim/solver/solver-base.hpp"

#include "lbm-sim/analysis/benchmarks.hpp"
#include "lbm-sim/analysis/error.hpp"
#include "lbm/logging.hpp"

#include "quill/LogMacros.h"

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
  void solve(SolverBase<dim, VelocitySet, cm_t, backend_t> &solver) {
    quill::Logger *simulation_logger =
        logging::create_or_get_logger("simulation");

    LOG_DEBUG(simulation_logger, "Initializing Simulation...");

    write_header(lattice.grid);

    solver.solve(lattice, params);

    LOG_DEBUG(simulation_logger, "Finished Simulation.");
  };

  /**
   * Errore rispetto a una soluzione analitica, sullo stile di
   * dealii::VectorTools: si passa la exact_solution (una Function<dim>,
   * es. CouetteSolution2D / PoiseuilleSolution2D, o una qualunque classe
   * derivata definita dall'utente) e il tipo di norma; il campo
   * approssimato (lattice.u) e la griglia sono presi internamente dalla
   * simulazione, non vanno passati dal chiamante.
   *
   * Equivalente di:
   *   VectorTools::integrate_difference(..., solution, exact_solution,
   *                                     error_per_cell, ..., norm_type);
   *   VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
   *
   * Ritorna l'errore globale assoluto (non normalizzato) nella norma
   * richiesta. Per l'errore relativo rispetto alla soluzione esatta,
   * vedi analysis::compute_error() in analysis/error.hpp.
   */
  double compute_error(const analysis::NormType &norm_type,
                       const analysis::Function<dim> &exact_solution) const {
    const auto error_per_cell =
        analysis::ErrorEvaluator<dim>::integrate_difference(
            lattice.grid, lattice.u, exact_solution);

    const double error = analysis::ErrorEvaluator<dim>::compute_global_error(
        error_per_cell, norm_type);

    return error;
  }

  /**
   * Errore rispetto al benchmark di Ghia et al. (1982), solo per la lid
   * cavity 2D: confronta lattice.u lungo le due centerline con le tabelle
   * di riferimento. lid_velocity e' la velocita' della parete mobile
   * (params.init_vel.dx), usata da Ghia per normalizzare i suoi dati.
   * norm_type sceglie la norma per ridurre i 17 punti tabulati a uno
   * scalare (default L2), stessa semantica di compute_error().
   *
   * Disponibile solo per dim == 2 -- la lid cavity di Ghia non ha un
   * equivalente 3D tabulato in questa libreria.
   */
  analysis::NormErrorResult compute_ghia_error(
      const std::string &filepath_in,
      const analysis::NormType norm_type = analysis::NormType::L2) const {
    if constexpr (dim == 2) {
      return analysis::compute_ghia_error(filepath_in, lattice, params.reyn_num,
                                          params.init_vel.dx);
    } else {
      static_assert(assertion::always_false<dim>,
                    "compute_ghia_error() is only defined for dim == 2");
    }
  }

  void output(const char *filepath,
              std::function<std::vector<double>(const Lattice<dim> &)>
                  extract_profile) {
    using namespace std::filesystem;

    quill::Logger *data_logger = logging::create_or_get_logger("data_log");

    path outpath(filepath);
    path parent = outpath.parent_path();

    if (!exists(parent)) {
      create_directories(parent);
    }

    std::ofstream fout(outpath, std::ios::binary);

    // FIXME: Should throw error?
    if (!fout.is_open()) {
      LOG_ERROR(data_logger, "Failed to create file: {}", filepath);
      return;
    }

    std::vector<double> profile = extract_profile(lattice);

    LOG_DEBUG(data_logger, "Writing header to file {}", filepath);

    std::string header = "%%profile " + collision_model_to_string(cm_t) + " " +
                         std::to_string(profile.size()) + " " +
                         std::to_string(params.init_vel.dx) + "\n";

    fout.write(header.data(), header.size());

    fout.write(reinterpret_cast<const char *>(profile.data()),
               profile.size() * sizeof(double));

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

    LOG_DEBUG(logging::create_or_get_logger("data_log"),
              "Writing norms header...");

    this->notifyListeners(std::move(buf));
  }
};

} // namespace lbm

#endif // __LBM_SIM_LBM_SIMULATION_HPP
