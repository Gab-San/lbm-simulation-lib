#pragma once

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/metadata.hpp"
#include "lbm-sim/types/base.hpp"

#include <sstream>

namespace lbm::format {

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

template <typename T, types::dim_t dim>
inline std::string file_format(utils::Point<T, dim> p) {
  std::ostringstream oss;
  if constexpr (dim == 2) {
    oss << p.x << "_" << p.y;
  } else {
    oss << p.x << "_" << p.y << "_" << p.z;
  }
  return oss.str();
}

inline std::string file_format(enum CollisionModel cm) {
  switch (cm) {
  case BGK:
    return "bgk";
  case TRT:
    return "trt";
  case MRT:
    return "mrt";
  }
  return std::to_string(cm);
}

inline std::string file_format(double reyn) {
  std::ostringstream oss;
  oss << reyn;
  return oss.str();
}

} // namespace lbm::format
