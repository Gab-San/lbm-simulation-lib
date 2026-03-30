#ifndef CONFIG_PARSER_HPP
#define CONFIG_PARSER_HPP

#include <vector>
#include <string>

#include "../defs.hpp"

namespace lbm_lbm {

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
struct Job {
	/// Number of cells along the x axis
	unsigned int grid_num_cells_x;

	/// Number of cells along the y axis
	unsigned int grid_num_cells_y;

	/// Number of iteration steps
	unsigned int num_steps;

	/// Number of frames
	///
	/// Frames contain the information about
	/// the norm of the velocity at a step t.
	unsigned int num_frames;

	/// Reynold number
	double reyn_num;

	/// Initial velocity of the lid
	double u_lid;

	/// Norms output file path 
	std::string vel_norm_out;

	/// Bench data output file path
	std::string vel_bench_out;

	/// Collision operator
	///
	/// Options:
	/// - BGK
	///
	/// Not supported:
	/// - TRT
	CollisionOperator coll_op;
};

/**
 * This function parses a configuration file
 * and tries to extract jobs from it.
 *
 * @param file_path path to config file
 */
std::vector<Job> parse_file(const std::string& file_path);

}

#endif
