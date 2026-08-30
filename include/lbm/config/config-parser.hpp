/**
 * @file config-parser.hpp
 * @brief Reading a TOML file into SimulationConfig objects.
 *
 * A file holds a top-level `conf` array, one entry per run, so a single file
 * can describe a sweep. Everything is validated as it is read and a failure
 * raises ConfigError naming the offending key -- there is no half-filled
 * config and no silently defaulted mandatory field.
 *
 * @see the "Configuration files" page for the schema and the error messages.
 */

#pragma once

#include "lbm/config/simulation-config.hpp"

#include <filesystem>
#include <iostream>
#include <vector>

namespace lbm::config {

template <unsigned short int dim>
std::vector<SimulationConfig<dim>>
parse_config(const std::filesystem::path &path);

inline void print_usage(const char *exe) {
  std::cerr << "Usage: " << exe << " <config.toml>\n\n"
            << "  config.toml   simulation configuration" << std::endl;
}

} // namespace lbm::config
