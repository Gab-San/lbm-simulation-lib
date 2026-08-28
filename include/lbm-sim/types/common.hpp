#pragma once

// The aliases below name Point/Vector as complete types, so pull in the
// definitions, not just types/fwd.hpp: every user of Coordinate<dim> et al.
// would otherwise have to remember to include core/point.hpp itself. Neither
// header includes this one, so there is no cycle.
#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"
#include "lbm-sim/types/base.hpp"
#include "lbm-sim/types/fwd.hpp"

#include <cstddef>

namespace lbm::types {

template <dim_t dim> using ValuePoint = utils::Point<double, dim>;

template <dim_t dim> using DimPoint = utils::Point<std::size_t, dim>;

template <dim_t dim> using Coordinate = utils::Point<int, dim>;

template <dim_t dim> using Direction = utils::Vector<int, dim>;

} // namespace lbm::types
