#include "ConfigParser.hpp"
#include <stdexcept>
#include <fstream>
#include <string>

#include <iostream>

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    auto end   = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

std::vector<std::string> split(std::string s, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;

    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
	token = s.substr(0, pos);
	tokens.push_back(token);
	s.erase(0, pos + delimiter.length());
    }

    tokens.push_back(s);
    return tokens;
}

std::string get_next_line(std::ifstream& in) {
    std::string line;
    if (!std::getline(in, line))
	throw std::runtime_error("[Error while parsing] Unexpected end of file");
    return trim(line);
}

using namespace lbm_lbm;

std::vector<Job> lbm_lbm::parse_file(const std::string& file_path) {

    std::vector<Job> jobs;
    
    std::ifstream in(file_path);
    
    if(!in)
	throw std::runtime_error("Cannot open file: " + file_path);

    std::string line;
    std::vector<std::string> toks;

    int parsing_state = 0;

    while(std::getline(in, line)) {
	// Skip all empty lines
	if(line.empty()) continue;
	line = trim(line);
	if(line[0] == '#') continue;

	Job currJob;
	while(parsing_state != 8) {
	
		toks = split(line, " ");
		switch(parsing_state) {

		    case 0:
			if(toks[0] != "Nx")
			    throw std::runtime_error("[Error while parsing] Nx line expected!");
			currJob.grid_num_cells_x = std::stoi(toks[1]);
			parsing_state = 1;
		    break; 

		    case 1:
			if(toks[0] != "Ny")
			    throw std::runtime_error("[Error while parsing] Ny line expected!");

			currJob.grid_num_cells_y = std::stoi(toks[1]);
			parsing_state = 2;
		    break;

		    case 2:
			if(toks[0] != "Re")
			    throw std::runtime_error("[Error while parsing] Re line expected!");
			
			currJob.reyn_num = std::stod(toks[1]);
			parsing_state = 3;
		    break;

		    case 3:
			if(toks[0] != "U_lid")
			    throw std::runtime_error("[Error while parsing] U_lid line expected!");
			
			currJob.u_lid = std::stod(toks[1]);
			parsing_state = 4;
		    break;

		    case 4:
			if(toks[0] != "Nsteps")
			    throw std::runtime_error("[Error while parsing] Nsteps line expected!");
			
			currJob.num_steps = std::stod(toks[1]);
			parsing_state = 5;
		    break;

		    case 5:
			if(toks[0] != "Nframes")
			    throw std::runtime_error("[Error while parsing] Nframes line expected!");
			
			currJob.num_frames = std::stod(toks[1]);
			parsing_state = 6;
		    break;

		    case 6:
			if(toks[0] != "out_norm")
			    throw std::runtime_error("[Error while parsing] out_norm line expected!");
			
			currJob.vel_norm_out = std::string(trim(toks[1]));
			parsing_state = 7;
		    break;

		    case 7:
			std::cout << currJob.vel_norm_out << std::endl;
			if(toks[0] != "out_bench")
			    throw std::runtime_error("[Error while parsing] out_bench line expected!");
			
			currJob.vel_bench_out = std::string(toks[1]);
			parsing_state = 8;
		   continue; 

		    default:
			throw std::runtime_error("[Error while parsing] Reached an undefined parsing state!");
		}

		do {

		    line = get_next_line(in);
		    line = trim(line);

		} while(line.empty() || line[0] == '#');
	}

	jobs.push_back(std::move(currJob));
	parsing_state = 0;
    }


    return jobs;
}
