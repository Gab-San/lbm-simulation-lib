#pragma once

#include "lbm-sim/types/base.hpp"
#include "lbm-sim/types/fwd.hpp"

#include <cstddef>

namespace lbm::types {

template <dim_t dim> using ValuePoint = utils::Point<double, dim>;

template <dim_t dim> using DimPoint = utils::Point<std::size_t, dim>;

template <dim_t dim> using Coordinate = utils::Point<int, dim>;

template <dim_t dim> using Direction = utils::Vector<int, dim>;

} // namespace lbm::types
