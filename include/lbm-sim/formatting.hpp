#pragma once

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/types/common.hpp"
#include <sstream>

namespace lbm::format {

inline std::string format_reyn(double reyn) {
  std::ostringstream oss;
  oss << reyn;
  return oss.str();
}

template <typename T, types::dim_t dim>
inline std::string csv_format(utils::Point<T, dim> p) {
  std::ostringstream oss;
  if constexpr (dim == 2) {
    oss << p.x << "x" << p.y;
  } else {
    oss << p.x << "x" << p.y << "x" << p.z;
  }
  return oss.str();
}

} // namespace lbm::format
