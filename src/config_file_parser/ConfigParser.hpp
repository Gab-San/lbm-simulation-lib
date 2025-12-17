#ifndef CONFIG_PARSER_HPP
#define CONFIG_PARSER_HPP

#include <vector>
#include <string>

namespace lbm_lbm {

struct Job {
	unsigned int grid_num_cells_x, grid_num_cells_y;
	unsigned int num_steps, num_frames;
	double reyn_num, u_lid;
	std::string vel_norm_out, vel_bench_out;
};

std::vector<Job> parse_file(const std::string& file_path);

}

#endif
