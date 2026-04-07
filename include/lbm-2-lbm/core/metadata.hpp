#pragma once

#include "lbm-2-lbm/types/core.hpp"

namespace lbm {
template<ExecutionBackend backend>
struct ExecutionContext;

template<>
struct ExecutionContext<ExecutionBackend::MPI> {
    int rank = 0;
    int world_size = 1;
};

template<>
struct ExecutionContext<ExecutionBackend::CUDA> {
    // Opaque stream handle to avoid forcing CUDA headers in common code.
    void* stream = nullptr;
    int block_x = 16;
    int block_y = 16;
    int block_z = 1;
};

template<int dim, enum types::CollisionModel coll_op>
struct Params;

template<>
struct Params<2, types::CollisionModel::BGK>  {
    const double reyn_num;
    const int num_cells_y;
    Params(
	const double _reyn_num,
	const double _num_cells_y
    ) : reyn_num(_reyn_num), num_cells_y(_num_cells_y) {};
};

template<>
struct Params<2, types::CollisionModel::TRT> {

};
}
