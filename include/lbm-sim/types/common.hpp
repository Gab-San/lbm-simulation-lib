#ifndef __LBM_SIM_TYPES_COMMON_HPP
#define __LBM_SIM_TYPES_COMMON_HPP

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/types/base.hpp"

namespace lbm {

namespace types {

template <types::dim_t dim> using ValuePoint = utils::Point<double, dim>;

template <types::dim_t dim> using DimPoint = utils::Point<std::size_t, dim>;

template <types::dim_t dim> using Coordinate = utils::Point<int, dim>;

} // namespace types

} // namespace lbm

#endif // __LBM_SIM_TYPES_COMMON_HPP
