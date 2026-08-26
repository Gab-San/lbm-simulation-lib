#pragma once
#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/types/common.hpp"

#include <string>

namespace lbm::config {

/// Parametri comuni a tutte le simulazioni, letti dal TOML.
///
/// NOTA: qui stanno solo gli scalari. La *geometria* (ostacoli, che parete
/// e' mobile, ecc.) resta nel main, perche' dipende dal tipo di problema:
/// e' il campo [problem].type del TOML a dire quale main deve prendersi in
/// carico un certo file di configurazione.
struct SimulationConfig {
  /// Nome della configurazione, ricavato dal nome del file .toml.
  /// Usato per comporre i path di output.
  std::string name;

  // [grid]
  unsigned int nx, ny;

  // [physics]
  CollisionModel collision;
  double reynolds;
  double init_vel_x;

  // [solver]
  unsigned int niters;
  unsigned int nframes;

  /// "openmp" | "cuda". Serve a distinguere i binari: senza, il main OpenMP
  /// e quello CUDA dello stesso problema matcherebbero gli stessi .toml e
  /// si ruberebbero le configurazioni a vicenda.
  std::string backend;

  // [output]
  /// Cartella di output (non un file): i nomi dei file vengono composti
  /// dal main a partire da `name`.
  std::string output_path;

  /// "<output_path>/norms_<name>.bin" - i frame delle norme di velocita'.
  std::string frames_file() const {
    return output_path + "/norms_" + name + ".bin";
  }

  /// "<output_path>/data_<name>.bin" - il profilo per il confronto coi
  /// benchmark.
  std::string profile_file() const {
    return output_path + "/data_" + name + ".bin";
  }
};

} // namespace lbm::config
