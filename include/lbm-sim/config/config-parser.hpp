#ifndef __LBM_SIM_CONFIG_CONFIG_PARSER_HPP
#define __LBM_SIM_CONFIG_CONFIG_PARSER_HPP

#include "lbm-sim/config/simulation-config.hpp"
#include "lbm-sim/metadata.hpp"

#include <toml++/toml.hpp>

#include <filesystem>
#include <iostream>
#include <string>

namespace lbm {
namespace config {

inline CollisionModel parse_collision_model(const std::string &s,
                                            const std::string &where) {
  if (s == "BGK") {
    return CollisionModel::BGK;
  }
  if (s == "TRT") {
    return CollisionModel::TRT;
  }
  if (s == "MRT") {
    return CollisionModel::MRT;
  }
  throw ConfigError(where + "operatore di collisione sconosciuto: '" + s +
                    "' (attesi BGK, TRT, MRT)");
}

inline ExecutionBackend parse_execution_backend(const std::string &s,
                                                const std::string &where) {
  if (s == "openmp") {
    return ExecutionBackend::OPEN_MP;
  }
  if (s == "cuda") {
    return ExecutionBackend::CUDA;
  }
  throw ConfigError(where + "Unknown backend: '" + s +
                    "' (wants openmp, cuda)");
}

// [solver].backend e' un enum, non una stringa: serve per i messaggi di
// errore di ensure_compatible. Se esiste gia' una funzione equivalente in
// backend.hpp usa quella al posto di questa.
inline std::string execution_backend_to_string(ExecutionBackend b) {
  switch (b) {
  case ExecutionBackend::OPEN_MP:
    return "openmp";
  case ExecutionBackend::CUDA:
    return "cuda";
  }
  return "unknown";
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
template <types::dim_t dim>
inline void parse_table(SimulationConfig<dim> &cfg, const toml::table &tbl,
                        const std::string &name = "") {
  const std::string where =
      name.empty() ? std::string("config") : ("config '" + name + "'");

  cfg.name = name;

  cfg.collision = parse_collision_model(
      tbl["collision"].value_or<std::string>("BGK"), where);
  cfg.backend = parse_execution_backend(
      tbl["backend"].value_or<std::string>("openmp"), where);

  // --- [lattice] ------------------------------------------------------------

  if (auto arr = tbl["lattice"]["size"].as_array()) {
    if (arr->size() != dim) {
      throw ConfigError(
          "config '" + cfg.name + "': e' un problema " + std::to_string(dim) +
          "D (dedotto " + (dim == 3 ? "dalla presenza" : "dall'assenza") +
          " di [grid].nz), ma questo binario e' " + std::to_string(dim) + "D");
    }
    std::size_t i = 0;
    for (auto &elem : *arr) {
      auto val = elem.template value<uint64_t>();
      if (!val)
        throw ConfigError(where + ": [lattice].size[" + std::to_string(i) +
                          "] mancante o non e' un intero");
      cfg.grid_size[i] = *val;
      ++i;
    }
  } else {
    throw ConfigError(where + ": [lattice].size e' obbligatorio");
  }

  // --- [physics] ---------------------------------------------------------

  {
    auto reynolds = tbl["physics"]["reynolds"].value<double>();
    if (!reynolds)
      throw ConfigError(where + ": [physics].reynolds e' obbligatorio");
    cfg.reynolds = *reynolds;
  }

  if (cfg.reynolds <= 0.0)
    throw ConfigError(where + ": [physics].reynolds deve essere > 0");

  std::array<double, dim> utmp;

  if (auto arr = tbl["physics"]["size"].as_array()) {
    if (arr->size() != dim) {
      throw ConfigError(
          "config '" + cfg.name + "': e' un problema " + std::to_string(dim) +
          "D (dedotto " + (dim == 3 ? "dalla presenza" : "dall'assenza") +
          " di [grid].nz), ma questo binario e' " + std::to_string(dim) + "D");
    }

    std::size_t i = 0;
    for (auto &elem : *arr) {
      auto val = elem.template value<double>();
      if (!val)
        throw ConfigError(where + ": [physics].size[" + std::to_string(i) +
                          "] mancante o non e' un numero");
      utmp[i] = *val;
      ++i;
    }
  } else {
    throw ConfigError(where + ": [physics].size e' obbligatorio");
  }

  cfg.u0 = utmp;

  // init_vel_x e' anche la velocita' caratteristica con cui
  // CollisionParams calcola nu = u*Ny/Re: se fosse 0 si avrebbe nu = 0.
  if (cfg.u0.dx <= 0.0)
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

  // --- [output] ----------------------------------------------------------

  {
    auto frames_out = tbl["output"]["frames"].value<std::string>();
    if (!frames_out)
      throw ConfigError(where + ": [output].frames e' obbligatorio");
    cfg.frames_out = *frames_out;

    auto profile_out = tbl["output"]["profile"].value<std::string>();
    if (!profile_out)
      throw ConfigError(where + ": [output].profile e' obbligatorio");
    cfg.profile_out = *profile_out;
  }

  // --- [backend] ---------------------------------------------------------
  // TODO: Add backend specifics
}

template <types::dim_t dim>
std::vector<SimulationConfig<dim>>
parse_config(const std::filesystem::path &path) {
  if (!std::filesystem::exists(path))
    throw ConfigError("file di configurazione inesistente: " + path.string());

  toml::table tbl;
  try {
    tbl = toml::parse_file(path.string());
  } catch (const toml::parse_error &err) {
    throw ConfigError("TOML malformato in " + path.string() + ": " +
                      std::string(err.description()));
  }

  const auto *arr = tbl["conf"].as_array();
  if (!arr)
    throw ConfigError("campo 'conf' mancante o non è un array in " +
                      path.string());

  std::vector<SimulationConfig<dim>> configs;
  configs.reserve(arr->size());

  for (const auto &elem : *arr) {
    const auto *sim_tbl = elem.as_table();
    if (!sim_tbl)
      throw ConfigError("elemento non valido in 'conf' in " + path.string());

    SimulationConfig<dim> cfg;
    parse_table(cfg, *sim_tbl, path.stem().string());
    configs.push_back(std::move(cfg));
  }

  return configs;
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
template <types::dim_t dim>
void ensure_compatible(const SimulationConfig<dim> &cfg,
                       CollisionModel collision,
                       const ExecutionBackend backend) {

  if (cfg.collision != collision)
    throw ConfigError("config '" + cfg.name + "': [physics].collision = '" +
                      collision_model_to_string(cfg.collision) +
                      "', ma questo binario e' compilato per '" +
                      collision_model_to_string(collision) + "'");

  if (cfg.backend != backend)
    throw ConfigError("config '" + cfg.name + "': [solver].backend = '" +
                      execution_backend_to_string(cfg.backend) +
                      "', ma questo binario e' '" +
                      execution_backend_to_string(backend) + "'");
}

inline void print_usage(const char *exe) {
  std::cerr << "Uso: " << exe << " <config.toml>\n\n"
            << "  config.toml   configurazione della simulazione" << std::endl;
}

} // namespace config
} // namespace lbm

#endif // __LBM_SIM_CONFIG_CONFIG_PARSER_HPP
