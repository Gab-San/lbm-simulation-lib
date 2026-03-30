#pragma once

#include "lbm-2-lbm/core/velocity_sets.hpp"
#include "lbm-2-lbm/core/collision_operators.hpp"

namespace lbm {
namespace types {
using BGK2D = BGKCollisionStrategy<2, D2Q9>;
}
}
