#ifndef __TYPES_COMMON_HPP
#define __TYPES_COMMON_HPP

#include "lbm-sim/core/point.hpp"

namespace lbm {
namespace types {
template <unsigned short int dim> using ValuePoint = utils::Point<double, dim>;

template <unsigned short int dim>
using DimPoint = utils::Point<std::size_t, dim>;

template <unsigned short int dim> using Coordinate = utils::Point<int, dim>;

} // namespace types
} // namespace lbm

#endif // __TYPES_COMMON_HPP
