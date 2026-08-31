#include "lbm/config/config-parser.hpp"

#include <toml++/toml.hpp>

#include <string>
namespace lbm::config {
template <unsigned short int dim>
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
          "D problem (deduced from the " + (dim == 3 ? "presence" : "absence") +
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

  if (auto arr = tbl["physics"]["init_vel"].as_array()) {
    if (arr->size() != dim) {
      throw ConfigError(
          "config '" + cfg.name + "': this is a " + std::to_string(dim) +
          "D problem (deduced from the " + (dim == 3 ? "presence" : "absence") +
          " of [grid].nz), but this binary is " + std::to_string(dim) + "D");
    }

    std::size_t i = 0;
    for (auto &elem : *arr) {
      auto val = elem.template value<double>();
      if (!val)
        throw ConfigError(where + ": [physics].size[" + std::to_string(i) +
                          "] is missing or is not a number");
      cfg.u0[i] = *val;
      ++i;
    }
  } else {
    throw ConfigError(where + ": [physics].size is mandatory");
  }

  // --- [solver] ----------------------------------------------------------
  cfg.niters = tbl["solver"]["niters"].value_or<unsigned int>(0);
  cfg.nframes = tbl["solver"]["nframes"].value_or<unsigned int>(0);

  if (cfg.niters == 0)
    throw ConfigError(where + ": [solver].niters is mandatory and must be > 0");

  if (cfg.nframes > cfg.niters)
    throw ConfigError(
        where + ": [solver].nframes (" + std::to_string(cfg.nframes) +
        ") cannot exceed [solver].niters (" + std::to_string(cfg.niters) + ")");

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

template <unsigned short int dim>
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

// Explicit template instantiations to satisfy the linker
template std::vector<SimulationConfig<2>>
parse_config<2>(const std::filesystem::path &);
template std::vector<SimulationConfig<3>>
parse_config<3>(const std::filesystem::path &);

} // namespace lbm::config
