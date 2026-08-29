#pragma once

#include <string>

namespace lbm {

enum ExecutionBackend { CUDA, OPEN_MP };

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

inline std::string backend_to_string(enum ExecutionBackend backend_t) {
  switch (backend_t) {
  case ExecutionBackend::OPEN_MP:
    return "OpenMP";
  case ExecutionBackend::CUDA:
    return "CUDA";
  }
  return std::to_string(backend_t);
}

} // namespace lbm
