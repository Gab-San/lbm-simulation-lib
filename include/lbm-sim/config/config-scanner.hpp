// config/config-scanner.hpp
#pragma once

#include "lbm-sim/config/config-parser.hpp"

#include <filesystem>
#include <functional>
#include <iostream>

namespace lbm::config {

/// Scansiona `configs_dir` ed esegue `on_match` solo sui file .toml il cui
/// campo [problem].type vale `wanted_type`.
inline void run_matching_configs(
    const std::filesystem::path &configs_dir, const std::string &wanted_type,
    const std::function<void(const SimulationConfig &, const toml::table &)>
        &on_match) {
  if (!std::filesystem::is_directory(configs_dir)) {
    std::cerr << "cartella di configurazione inesistente: " << configs_dir
              << "\n";
    return;
  }

  for (const auto &entry : std::filesystem::directory_iterator(configs_dir)) {
    if (entry.path().extension() != ".toml")
      continue;

    toml::table tbl;
    try {
      tbl = toml::parse_file(entry.path().string());
    } catch (const toml::parse_error &err) {
      std::cerr << "config malformato, scartato: " << entry.path() << " ("
                << err.description() << ")\n";
      continue;
    }

    const std::string type = tbl["problem"]["type"].value_or<std::string>("");
    if (type != wanted_type) {
      continue; // non e' per questo main, scartato silenziosamente
    }

    try {
      on_match(parse_common(tbl, entry.path().stem().string()), tbl);
    } catch (const std::exception &err) {
      std::cerr << "config non valido, scartato: " << entry.path() << " ("
                << err.what() << ")\n";
    }
  }
}

} // namespace lbm::config
