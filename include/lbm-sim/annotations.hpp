/**
 * @file annotations.hpp
 * @brief The portable spelling of the @c restrict qualifier.
 *
 * @c restrict is C, not C++, and every compiler exposes it under a different
 * name. The macro is used on the population pointers handed to the collision
 * kernels, where promising the compiler that the buffers do not alias is
 * what allows the inner loop to vectorise.
 *
 * Expands to nothing on a compiler that offers no equivalent, so the code
 * stays correct and merely loses the optimisation.
 *
 * @see backend/cuda/annotations.hpp for @c LBM_HD_FUNC and
 *      backend/omp/annotations.hpp for @c UNROLL_FULL.
 */

#pragma once

#ifndef RESTRICT
#if defined(_MSC_VER)
#define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define RESTRICT restrict
#else
#define RESTRICT
#endif
#endif
