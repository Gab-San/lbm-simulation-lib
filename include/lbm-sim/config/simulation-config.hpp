/**
 * @file simulation-config.hpp
 * @brief SimulationConfig: one run's worth of settings, and the error type
 *        that reports a bad one.
 *
 * The struct is a plain data carrier -- the parsing and the validation live
 * in config-parser.hpp. It holds only what the file actually supplies: the
 * dimension, the collision operator and the backend are compile-time
 * parameters of the binary, not fields.
 *
 * @see the "Configuration files" page for the TOML schema.
 */

#pragma once

#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/base.hpp"
#include <array>
#include <stdexcept>
#include <string>

namespace lbm::config {

/**
 * \brief Configuration error.
 *
 * Covers everything that makes a configuration unusable *before* the
 * simulation starts: unreadable file, malformed TOML, missing mandatory
 * field, out-of-range value, or a configuration that is self-consistent but
 * incompatible with the binary reading it (dimension, problem type,
 * collision operator, backend).
 *
 * It is a dedicated type rather than a plain std::runtime_error precisely so
 * that the mains can tell "the user passed the wrong config" (print the
 * message and exit with 1) apart from an error raised during the solve.
 */
class ConfigError : public std::runtime_error {
public:
  explicit ConfigError(const std::string &what) : std::runtime_error(what) {}
};

/**
 * \brief Configuration of *one* simulation, one entry of the `conf` array of
 * a .toml file.
 *
 * The problem dimension is not a field: parse_config<dim>() is instantiated
 * with the dimension the binary was compiled for, and an array whose length
 * does not match is rejected. See config-parser.hpp and the "Configuration
 * files" page.
 *
 * @verbatim
   [[conf]]

   [conf.lattice]
   size = [129, 129]        # [nx, ny] in 2D, [nx, ny, nz] in 3D

   [conf.physics]
   reynolds = 100.0
   size     = [0.1, 0.0]    # reference velocity, size[0] must be > 0

   [conf.solver]
   niters  = 10000
   nframes = 100

   [conf.output]
   frames  = "out/frames.bin"
   profile = "out/profile.dat"

   [conf.backend]
   n_threads = 8            # optional
   @endverbatim
 */
template <types::dim_t dim> struct SimulationConfig {

  /// Grid cells. `nz` is 1 (and must not be read) when dim == 2.
  std::array<uint64_t, dim> grid_size;

  /// Stem of the configuration file name, filled in by the parser and used
  /// in error messages and output basenames.
  std::string name;
  /// Reynolds number.
  double reynolds = 100.0;

  /// Reference velocity. `init_vel_z` exists only for dim == 3.
  ///
  /// \note init_vel_x is not just the initial velocity: it is the
  /// characteristic velocity CollisionParams uses to compute nu = u*Ny/Re,
  /// so it must be strictly positive even in problems (Poiseuille) where no
  /// wall moves.
  utils::Vector<double, dim> u0;

  /// Number of solver iterations.
  unsigned int niters = 0;

  /// Number of frames saved during the solve. Each frame holds the velocity
  /// norms at a given time step.
  unsigned int nframes = 0;

  /// Output file for the norms, from [conf.output].frames.
  std::string frames_out;

  /// Output file for the velocity profile, from [conf.output].profile.
  std::string profile_out;

  /// Thread count requested by [conf.backend].n_threads; 0 means "let the
  /// backend decide".
  unsigned int n_threads;
};

} // namespace lbm::config
