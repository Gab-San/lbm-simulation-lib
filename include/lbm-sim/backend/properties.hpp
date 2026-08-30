/**
 * @file properties.hpp
 * @brief Umbrella header for the BackendProperties specialisations.
 */

#pragma once

#include "lbm-sim/metadata.hpp"

/// Umbrella header: pulls in every `BackendProperties` specialization that is
/// available in this translation unit. Include this to *configure* a backend.
///
/// The axis to switch on is what the build enabled, not which compiler is
/// parsing the file. A `.cu` translation unit still runs plenty of work on the
/// host -- computing the boundary mask, for instance -- so it wants the OpenMP
/// specialization just as much as a `.cpp` one does. Hence `_OPENMP` (defined
/// by the host compiler whenever OpenMP is on) rather than
/// `!defined(__CUDACC__)`: under nvcc both branches can, and should, be taken.
///
/// Code that only needs to *name* the template should include
/// "lbm-sim/backend/fwd.hpp" instead and stay free of <omp.h> and the CUDA
/// runtime.

#if defined(_OPENMP)
#include "lbm-sim/backend/omp/properties.hpp"
#endif

// __CUDACC__ guards *parsing*: only nvcc can read the .cuh. If the CUDA
// properties end up needing nothing but the host-side runtime API
// (cudaSetDevice, cudaGetDeviceProperties, ...), rename it to .hpp and switch
// this to LBM_ENABLE_CUDA -- a plain .cpp could then configure the GPU too.
// Note that LBM_ENABLE_CUDA is currently a CMake option only: the
// target_compile_definitions that would make it visible to the preprocessor is
// commented out in cmake/Cuda.cmake.
#if defined(__CUDACC__)
#include "lbm-sim/backend/cuda/properties.cuh"
#endif
