/**
 * @file fwd.hpp
 * @brief Forward declarations of Point and Vector.
 *
 * Lets a header name @c Point<T,dim> or @c Vector<T,dim> in a signature
 * without pulling in their definitions -- and, more importantly, lets
 * point.hpp and vector.hpp refer to each other (a Vector is built from two
 * Points) without a circular include.
 *
 * Include types/common.hpp instead if the complete types are needed.
 */

#pragma once

#include "lbm-sim/types/base.hpp"

namespace lbm::utils {

template <typename T, types::dim_t dim> struct Point;
template <typename T, types::dim_t dim> struct Vector;

} // namespace lbm::utils
