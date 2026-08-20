#ifndef __TYPES_COMMON_HPP
#define __TYPES_COMMON_HPP

#include "lbm-sim/core/point.hpp"

#include <cstdint>
#include <vector>

namespace lbm {
namespace types {
template <unsigned short int dim> using ValuePoint = utils::Point<double, dim>;

template <unsigned short int dim>
using DimPoint = utils::Point<std::size_t, dim>;

template <unsigned short int dim> using Coordinate = utils::Point<int, dim>;

using boundary_t = uint8_t;
using boundary_mask_t = std::vector<boundary_t>;

} // namespace types
} // namespace lbm

#endif // __TYPES_COMMON_HPP
