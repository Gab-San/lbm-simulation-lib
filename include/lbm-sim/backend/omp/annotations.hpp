/**
 * @file annotations.hpp
 * @brief @c UNROLL_FULL: the portable spelling of a full-unroll pragma.
 *
 * Applied to the loops over the @c ndir directions of a velocity set. Their
 * trip count is a compile-time constant (9, 19 or 27) and the body indexes
 * @c constexpr tables, so unrolling turns the table lookups into immediates.
 *
 * Expands to nothing on a compiler with no equivalent pragma, MSVC included.
 */

#ifndef __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP
#define __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP

#pragma once
#if defined(__clang__)
#define UNROLL_FULL _Pragma("clang loop unroll(full)")
#elif defined(__GNUC__)
#define UNROLL_FULL _Pragma("GCC unroll 65534")
#else
#define UNROLL_FULL
#endif

#endif // __LBM_SIM_BACKEND_OMP_ANNOTATIONS_HPP
