#pragma once

#include <string>

namespace lbm {

// FIXME: DECOMMENT AND REWIRE
// enum ExecutionBackend { CUDA, OPEN_MP };

enum CollisionModel { BGK, TRT, MRT };

inline std::string collision_model_to_string(enum CollisionModel cm_t) {
  switch (cm_t) {
  case CollisionModel::BGK:
    return "BGK";
  case CollisionModel::TRT:
    return "TRT";
  case CollisionModel::MRT:
    return "MRT";
  }
  return std::to_string(cm_t);
}

} // namespace lbm
