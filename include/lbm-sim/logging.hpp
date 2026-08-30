/**
 * @file logging.hpp
 * @brief The logging entry point for lbm-sim.
 *
 * Forwards to the backend-agnostic facade in lbm/logging/logging.hpp, and
 * additionally registers the formatters for Point and Vector when the quill
 * backend is in use, so `LBM_LOG_INFO(log, "grid {}", grid.size)` works
 * without the call site knowing which backend it compiled against.
 */

#pragma once
#include "lbm/logging/logging.hpp"

#if LBM_LOG_BACKEND == LBM_LOG_BACKEND_QUILL
#include "lbm-sim/logging/formatters.hpp"
#endif
