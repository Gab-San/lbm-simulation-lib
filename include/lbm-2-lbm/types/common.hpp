#pragma once

#include "lbm-2-lbm/core/point.hpp"
#include <cstddef>

namespace lbm {
namespace types {
template<int Pdim>
using ValuePoint = utils::Point<double, Pdim>;

template<int Pdim>
using IndexPoint = utils::Point<unsigned int, Pdim>;

template<int Pdim>
using DimPoint = utils::Point<std::size_t, Pdim>;
}
}

