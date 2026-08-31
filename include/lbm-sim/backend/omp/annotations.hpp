#ifndef __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP
#define __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP

// Manage loop unrolling pragmas for different compilers
#pragma once
#if defined(__clang__)
#define UNROLL_FULL _Pragma("clang loop unroll(full)")
#elif defined(__GNUC__)
#define UNROLL_FULL _Pragma("GCC unroll 65534")
#else
#define UNROLL_FULL
#endif

#endif // __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP
