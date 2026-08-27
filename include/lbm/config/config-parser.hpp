#ifndef __LBM_SIM_CONFIG_CONFIG_PARSER_HPP
#define __LBM_SIM_CONFIG_CONFIG_PARSER_HPP

#include "lbm-sim/config/simulation-config.hpp"

#include <toml++/toml.hpp>

#include <filesystem>
#include <string>

namespace lbm {
namespace config {

inline CollisionModel parse_collision_model(const std::string &s) {
  if (s == "BGK")
    return CollisionModel::BGK;
  if (s == "TRT")
    return CollisionModel::TRT;
  if (s == "MRT")
    return CollisionModel::MRT;
  throw ConfigError("operatore di collisione sconosciuto: '" + s +
                    "' (attesi BGK, TRT, MRT)");
}

/**
 * \brief Legge una SimulationConfig da una tabella TOML gia' parsata.
 *
 * I due path di output non stanno nel TOML: li passa il main dalla riga di
 * comando, perche' identificano *questa* esecuzione, non la configurazione
 * fisica (lo stesso .toml puo' essere rilanciato scrivendo altrove).
 *
 * La dimensione si deduce da [grid].nz: assente -> 2D, presente -> 3D.
 * [grid].dim e' opzionale e serve solo a far ricontrollare la deduzione.
 */
inline SimulationConfig parse_common(const toml::table &tbl,
                                     const std::string &name = "") {
  SimulationConfig cfg;
  cfg.name = name;

  const std::string where = name.empty() ? std::string("config")
                                         : ("config '" + name + "'");

  // --- [problem] ---------------------------------------------------------
  // Opzionale: ogni binario sa gia' quale problema risolve, quindi il tipo
  // nel file non serve a sceglierlo. Se pero' c'e', ensure_compatible() lo
  // verifica, cosi' passare per sbaglio la config di un altro problema
  // viene segnalato invece di girare in silenzio con la geometria sbagliata.
  cfg.problem = tbl["problem"]["type"].value_or<std::string>("");

  // --- [grid] ------------------------------------------------------------
  const auto nx = tbl["grid"]["nx"].value<int64_t>();
  const auto ny = tbl["grid"]["ny"].value<int64_t>();
  if (!nx || !ny)
    throw ConfigError(where + ": [grid].nx e [grid].ny sono obbligatori");
  if (*nx <= 0 || *ny <= 0)
    throw ConfigError(where + ": [grid].nx e [grid].ny devono essere > 0");

  cfg.nx = static_cast<unsigned int>(*nx);
  cfg.ny = static_cast<unsigned int>(*ny);

  // nz presente => problema 3D.
  const auto nz = tbl["grid"]["nz"].value<int64_t>();
  if (nz) {
    if (*nz <= 0)
      throw ConfigError(where + ": [grid].nz deve essere > 0");
    cfg.dim = 3;
    cfg.nz = static_cast<unsigned int>(*nz);
  } else {
    cfg.dim = 2;
    cfg.nz = 1;
  }

  // [grid].dim e' ridondante: se c'e', deve confermare la deduzione.
  if (const auto declared = tbl["grid"]["dim"].value<int64_t>()) {
    if (*declared != static_cast<int64_t>(cfg.dim))
      throw ConfigError(
          where + ": [grid].dim = " + std::to_string(*declared) +
          " ma dalla presenza/assenza di [grid].nz risulta un problema " +
          std::to_string(cfg.dim) + "D");
  }

  // --- [physics] ---------------------------------------------------------
  cfg.collision = parse_collision_model(
      tbl["physics"]["collision"].value_or<std::string>("BGK"));

  cfg.reynolds = tbl["physics"]["reynolds"].value_or(100.0);
  if (cfg.reynolds <= 0.0)
    throw ConfigError(where + ": [physics].reynolds deve essere > 0");

  cfg.init_vel_x = tbl["physics"]["init_vel_x"].value_or(0.0);
  cfg.init_vel_y = tbl["physics"]["init_vel_y"].value_or(0.0);
  cfg.init_vel_z = tbl["physics"]["init_vel_z"].value_or(0.0);

  // init_vel_x e' anche la velocita' caratteristica con cui
  // CollisionParams calcola nu = u*Ny/Re: se fosse 0 si avrebbe nu = 0.
  if (cfg.init_vel_x <= 0.0)
    throw ConfigError(where +
                      ": [physics].init_vel_x deve essere > 0 (e' la "
                      "velocita' caratteristica usata per nu = u*Ny/Re)");

  // --- [solver] ----------------------------------------------------------
  cfg.niters = tbl["solver"]["niters"].value_or<unsigned int>(0);
  cfg.nframes = tbl["solver"]["nframes"].value_or<unsigned int>(0);
  if (cfg.niters == 0)
    throw ConfigError(where + ": [solver].niters e' obbligatorio e > 0");
  if (cfg.nframes > cfg.niters)
    throw ConfigError(where + ": [solver].nframes (" +
                      std::to_string(cfg.nframes) +
                      ") non puo' superare [solver].niters (" +
                      std::to_string(cfg.niters) + ")");

  cfg.backend = tbl["solver"]["backend"].value_or<std::string>("openmp");

  return cfg;
}

/// Parsa un singolo file .toml. `name` viene ricavato dal nome del file.
inline SimulationConfig parse_config(const std::filesystem::path &path) {
  if (!std::filesystem::exists(path))
    throw ConfigError("file di configurazione inesistente: " + path.string());

  toml::table tbl;
  try {
    tbl = toml::parse_file(path.string());
  } catch (const toml::parse_error &err) {
    throw ConfigError("TOML malformato in " + path.string() + ": " +
                      std::string(err.description()));
  }

  return parse_common(tbl, path.stem().string());
}

/**
 * \brief Verifica che la configurazione sia 2D o 3D come si aspetta il
 * chiamante.
 *
 * La dimensione non e' dichiarata: la deduce parse_common() dalla presenza
 * di [grid].nz. Questo e' il controllo che impedisce di passare una config
 * 3D a un binario 2D (e viceversa), che altrimenti leggerebbe nx/ny e
 * ignorerebbe silenziosamente la terza dimensione.
 */
inline void ensure_dim(const SimulationConfig &cfg, types::dim_t dim) {
  if (cfg.dim != dim)
    throw ConfigError("config '" + cfg.name + "': e' un problema " +
                      std::to_string(cfg.dim) +
                      "D (dedotto " + (cfg.dim == 3 ? "dalla presenza"
                                                    : "dall'assenza") +
                      " di [grid].nz), ma questo binario e' " +
                      std::to_string(dim) + "D");
}

/**
 * \brief Verifica che la configurazione sia eseguibile da *questo* binario.
 *
 * Prima ogni main scandiva una cartella e scartava in silenzio i .toml non
 * suoi; adesso che il file lo sceglie l'utente sulla riga di comando, un
 * file non compatibile e' un errore da segnalare, non da ignorare.
 *
 * [problem].type e' opzionale: se il file non lo dichiara non si controlla
 * nulla, perche' il binario sa gia' quale problema risolve.
 */
inline void ensure_compatible(const SimulationConfig &cfg,
                              const std::string &problem,
                              CollisionModel collision,
                              const std::string &backend, types::dim_t dim) {
  if (!cfg.problem.empty() && cfg.problem != problem)
    throw ConfigError("config '" + cfg.name + "': [problem].type = '" +
                      cfg.problem + "', ma questo binario esegue '" + problem +
                      "'");

  if (cfg.collision != collision)
    throw ConfigError("config '" + cfg.name + "': [physics].collision = '" +
                      collision_model_to_string(cfg.collision) +
                      "', ma questo binario e' compilato per '" +
                      collision_model_to_string(collision) + "'");

  if (cfg.backend != backend)
    throw ConfigError("config '" + cfg.name + "': [solver].backend = '" +
                      cfg.backend + "', ma questo binario e' '" + backend +
                      "'");

  ensure_dim(cfg, dim);
}

} // namespace config
} // namespace lbm

#endif // __LBM_SIM_CONFIG_CONFIG_PARSER_HPP
