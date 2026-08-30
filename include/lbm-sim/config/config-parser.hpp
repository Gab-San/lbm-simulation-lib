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
  throw ConfigError(where + "unknown collision operator: '" + s +
                    "' (expected BGK, TRT, MRT)");
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

// [solver].backend is an enum, not a string: this is needed for the error
// messages of ensure_compatible. If an equivalent function already exists in
// backend.hpp, use that one instead.
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
 * \brief Reads a SimulationConfig from an already parsed TOML table.
 *
 * Reads [lattice].size, [physics].reynolds, [physics].size, [solver].niters,
 * [solver].nframes, [output].frames, [output].profile and
 * [backend].n_threads. Any other key in the table is ignored.
 *
 * The dimension is not read from the file: it is the template argument, and
 * an array whose length does not match it is rejected with a ConfigError
 * naming both.
 *
 * @param cfg  Destination, filled in place.
 * @param tbl  One entry of the `conf` array.
 * @param name Run name, used in the error messages.
 */
template <types::dim_t dim>
inline void parse_table(SimulationConfig<dim> &cfg, const toml::table &tbl,
                        const std::string &name = "") {
  const std::string where =
      name.empty() ? std::string("config") : ("config '" + name + "'");

  cfg.name = name;

  // --- [lattice] ------------------------------------------------------------

  if (auto arr = tbl["lattice"]["size"].as_array()) {
    if (arr->size() != dim) {
      throw ConfigError(
          "config '" + cfg.name + "': this is a " + std::to_string(dim) +
          "D problem (deduced from the " +
          (dim == 3 ? "presence" : "absence") +
          " of [grid].nz), but this binary is " + std::to_string(dim) + "D");
    }
    std::size_t i = 0;
    for (auto &elem : *arr) {
      auto val = elem.template value<uint64_t>();
      if (!val)
        throw ConfigError(where + ": [lattice].size[" + std::to_string(i) +
                          "] is missing or is not an integer");
      cfg.grid_size[i] = *val;
      ++i;
    }
  } else {
    throw ConfigError(where + ": [lattice].size is mandatory");
  }

  // --- [physics] ---------------------------------------------------------

  {
    auto reynolds = tbl["physics"]["reynolds"].value<double>();
    if (!reynolds)
      throw ConfigError(where + ": [physics].reynolds is mandatory");
    cfg.reynolds = *reynolds;
  }

  if (cfg.reynolds <= 0.0)
    throw ConfigError(where + ": [physics].reynolds must be > 0");

  std::array<double, dim> utmp;

  if (auto arr = tbl["physics"]["size"].as_array()) {
    if (arr->size() != dim) {
      throw ConfigError(
          "config '" + cfg.name + "': this is a " + std::to_string(dim) +
          "D problem (deduced from the " +
          (dim == 3 ? "presence" : "absence") +
          " of [grid].nz), but this binary is " + std::to_string(dim) + "D");
    }

    std::size_t i = 0;
    for (auto &elem : *arr) {
      auto val = elem.template value<double>();
      if (!val)
        throw ConfigError(where + ": [physics].size[" + std::to_string(i) +
                          "] is missing or is not a number");
      utmp[i] = *val;
      ++i;
    }
  } else {
    throw ConfigError(where + ": [physics].size is mandatory");
  }

  cfg.u0 = utmp;

  // init_vel_x is also the characteristic velocity CollisionParams uses to
  // compute nu = u*Ny/Re: were it 0, nu would be 0 as well.
  if (cfg.u0.dx <= 0.0)
    throw ConfigError(where +
                      ": [physics].init_vel_x must be > 0 (it is the "
                      "characteristic velocity used for nu = u*Ny/Re)");

  // --- [solver] ----------------------------------------------------------
  cfg.niters = tbl["solver"]["niters"].value_or<unsigned int>(0);
  cfg.nframes = tbl["solver"]["nframes"].value_or<unsigned int>(0);

  if (cfg.niters == 0)
    throw ConfigError(where + ": [solver].niters is mandatory and must be > 0");

  if (cfg.nframes > cfg.niters)
    throw ConfigError(where + ": [solver].nframes (" +
                      std::to_string(cfg.nframes) +
                      ") cannot exceed [solver].niters (" +
                      std::to_string(cfg.niters) + ")");

  // --- [output] ----------------------------------------------------------

  {
    auto frames_out = tbl["output"]["frames"].value<std::string>();
    if (!frames_out)
      throw ConfigError(where + ": [output].frames is mandatory");
    cfg.frames_out = *frames_out;

    auto profile_out = tbl["output"]["profile"].value<std::string>();
    if (!profile_out)
      throw ConfigError(where + ": [output].profile is mandatory");
    cfg.profile_out = *profile_out;
  }

  // --- [backend] ---------------------------------------------------------

  {
    cfg.n_threads = tbl["backend"]["n_threads"].value_or<unsigned int>(0);
  }
}

template <types::dim_t dim>
std::vector<SimulationConfig<dim>>
parse_config(const std::filesystem::path &path) {
  if (!std::filesystem::exists(path))
    throw ConfigError("configuration file does not exist: " + path.string());

  toml::table tbl;
  try {
    tbl = toml::parse_file(path.string());
  } catch (const toml::parse_error &err) {
    throw ConfigError("malformed TOML in " + path.string() + ": " +
                      std::string(err.description()));
  }

  const auto *arr = tbl["conf"].as_array();
  if (!arr)
    throw ConfigError("field 'conf' is missing or is not an array in " +
                      path.string());

  std::vector<SimulationConfig<dim>> configs;
  configs.reserve(arr->size());

  for (const auto &elem : *arr) {
    const auto *sim_tbl = elem.as_table();
    if (!sim_tbl)
      throw ConfigError("invalid element in 'conf' in " + path.string());

    SimulationConfig<dim> cfg;
    parse_table(cfg, *sim_tbl, path.stem().string());
    configs.push_back(std::move(cfg));
  }

  return configs;
}

/**
 * \brief Checks that the configuration can be run by *this* binary.
 *
 * Each main used to scan a directory and silently discard the .toml files
 * that were not its own; now that the user picks the file on the command
 * line, an incompatible file is an error to report, not to ignore.
 *
 * [problem].type is optional: if the file does not declare it nothing is
 * checked, because the binary already knows which problem it solves.
 */
template <types::dim_t dim>
void ensure_compatible(const SimulationConfig<dim> &cfg,
                       CollisionModel collision,
                       const ExecutionBackend backend) {

  if (cfg.collision != collision)
    throw ConfigError("config '" + cfg.name + "': [physics].collision = '" +
                      collision_model_to_string(cfg.collision) +
                      "', but this binary is compiled for '" +
                      collision_model_to_string(collision) + "'");

  if (cfg.backend != backend)
    throw ConfigError("config '" + cfg.name + "': [solver].backend = '" +
                      execution_backend_to_string(cfg.backend) +
                      "', but this binary is '" +
                      execution_backend_to_string(backend) + "'");
}

inline void print_usage(const char *exe) {
  std::cerr << "Usage: " << exe << " <config.toml>\n\n"
            << "  config.toml   simulation configuration" << std::endl;
}

} // namespace config
} // namespace lbm
