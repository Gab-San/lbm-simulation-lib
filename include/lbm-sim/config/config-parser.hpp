#pragma once
#include "lbm-sim/config/simulation-config.hpp"

#include <toml++/toml.hpp>

#include <filesystem>
#include <stdexcept>
#include <string>

namespace lbm::config {

inline CollisionModel parse_collision_model(const std::string &s) {
  if (s == "BGK")
    return CollisionModel::BGK;
  if (s == "TRT")
    return CollisionModel::TRT;
  if (s == "MRT")
    return CollisionModel::MRT;
  throw std::invalid_argument("collision model sconosciuto: " + s);
}

/// Legge i campi comuni da una tabella TOML gia' parsata.
///
inline SimulationConfig parse_common(const toml::table &tbl,
                                     const std::string &name = "") {
  SimulationConfig cfg;
  cfg.name = name;

  // value<T>() + controllo esplicito per i campi obbligatori;
  // value_or() per quelli che hanno un default sensato.
  auto nx = tbl["grid"]["nx"].value<int64_t>();
  auto ny = tbl["grid"]["ny"].value<int64_t>();
  if (!nx || !ny)
    throw std::runtime_error("config: [grid].nx/ny obbligatori");
  cfg.nx = static_cast<unsigned int>(*nx);
  cfg.ny = static_cast<unsigned int>(*ny);

  auto collision_str = tbl["physics"]["collision"].value<std::string>();
  cfg.collision = parse_collision_model(collision_str.value_or("BGK"));

  cfg.reynolds = tbl["physics"]["reynolds"].value_or(100.0);
  cfg.init_vel_x = tbl["physics"]["init_vel_x"].value_or(0.1);

  cfg.niters = tbl["solver"]["niters"].value_or<unsigned int>(1000);
  cfg.nframes = tbl["solver"]["nframes"].value_or<unsigned int>(100);
  cfg.backend = tbl["solver"]["backend"].value_or<std::string>("openmp");

  cfg.output_path = tbl["output"]["path"].value_or<std::string>("./out");

  return cfg;
}

/// Parsa un singolo file di configurazione. `name` viene ricavato dal nome
/// del file (senza estensione).
inline SimulationConfig parse_config(const std::filesystem::path &path) {
  toml::table tbl;
  try {
    tbl = toml::parse_file(path.string());
  } catch (const toml::parse_error &err) {
    throw std::runtime_error("errore di parsing TOML in " + path.string() +
                             ": " + std::string(err.description()));
  }

  return parse_common(tbl, path.stem().string());
}

} // namespace lbm::config
