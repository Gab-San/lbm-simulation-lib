/**
 * @file annotations.hpp
 * @brief @c LBM_HD_FUNC: the host/device qualifier, or nothing.
 *
 * Under nvcc it expands to @c __host__ @c __device__, so the function is
 * compiled for both sides; under any other compiler it disappears, and the
 * same header remains valid C++ in a plain @c .cpp translation unit.
 *
 * This is what lets one implementation of the grid indexing, the vector
 * algebra, the boundary conditions and the collision kernels serve both
 * backends. A function marked with it must therefore stay device-compatible:
 * no allocation, no exceptions, no standard containers, no iostreams.
 */

#ifndef __LBM_SIM_BACKEND_CUDA_METADATA_HPP
#define __LBM_SIM_BACKEND_CUDA_METADATA_HPP

#ifdef __CUDACC__
#define LBM_HD_FUNC __host__ __device__
#else
#define LBM_HD_FUNC
#endif

#endif // __LBM_SIM_BACKEND_CUDA_METADATA_HPP
