#ifndef CONFIG_PARSER_HPP
#define CONFIG_PARSER_HPP

#include "lbm-sim/collision-operators/metadata.hpp"
#include "lbm-sim/core/types.hpp"
#include "lbm-sim/core/vector.hpp"

#include <string>
#include <vector>

namespace lbm {
namespace config {

/**
 * \brief This struct represents a job.
 * A job is an instance of a simulation.
 *
 * The parameters of a lid cavity simulation are:
 * - grid_size (num_cells_x, num_cells_y)
 * - reyn_num: reynold number
 * - u_lid: initial velocity of the lid
 * - coll_op: collision operator
 * - num_steps: number of iteration steps
 * - num_frames: number of frames to save
 *
 *
 * Frames contain the information about
 * the norm of the velocity at a step t.
 *
 *
 * Jobs are divided from one another by
 * whitelines.
 * There is no special string of characters to divide jobs,
 * but all of the parameters of a job must be configured.
 *
 * \note TRT operator is not supported
 *
 *
 * \section Configuration File
 *
 * The configuration file should look as following:
 *
 * @verbatim
  # Line Comment
  Nx 100
  Ny 100
  Re 100
  U_lid 0.1
  Nsteps 5000
  Nframes 50
  out_norm out/norm.txt
  out_bench out/data.txt
  Col BGK
  @endverbatim
 *
 */
template <types::dim_t dim> struct Config {

  /// Size of the grid.
  types::DimPoint<dim> grid_size;

  /// Number of iteration steps
  const unsigned int iters;

  /// Number of frames
  ///
  /// Frames contain the information about
  /// the norm of the velocity at a step t.
  const unsigned int frames;

  /// Reynold number
  const double reyn_num;

  /// Initial velocity of the fluid
  const utils::Vector<double, 2> init_vel;

  /// Norms output file path
  const std::string frames_filepath_out;

  /// Bench data output file path
  const std::string profile_filepath_out;

  /// Collision operator
  ///
  /// Options:
  /// - BGK
  /// - TRT
  CollisionModel collision_operator;

  Config(const types::DimPoint<dim> grid_size_, const unsigned int iters_,
         const unsigned int frames_, const double reyn_num_,
         const utils::Vector<double, dim> init_vel_,
         const std::string frames_filepath_out_,
         const std::string profile_filepath_out_)
      : grid_size(grid_size_), iters(iters_), frames(frames_),
        reyn_num(reyn_num_), init_vel(init_vel_),
        frames_filepath_out(frames_filepath_out_),
        profile_filepath_out(profile_filepath_out_) {}
};

/**
 * This function parses a configuration file
 * and tries to extract jobs from it.
 *
 * @param file_path path to config file
 */
template <types::dim_t dim>
std::vector<Config<dim>> parse_file(const std::string &file_path);

} // namespace config
} // namespace lbm

#endif
