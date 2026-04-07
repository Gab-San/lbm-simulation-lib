#pragma once

namespace lbm {
enum ExecutionBackend {
    CUDA,
    MPI
};

namespace types {

enum CollisionModel {
    BGK,
    TRT,
    MRT
};

}
}