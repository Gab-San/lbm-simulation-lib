#pragma once

// Quill backend. Included by lbm/logging/logging.hpp; never include directly.
//
// Owes the facade exactly two things: the definition of `Logger`, and the
// eight LBM_LOG_* macros.

#include "lbm/logging/quill/nvcc-compat.hpp"

#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/Sink.h"

#include <memory>

namespace lbm::logging {

// A wrapper rather than `using Logger = quill::Logger`. Costs one pointer
// load per statement and buys the guarantee that core code cannot call a
// quill method on a logger and only find out when the backend is switched.
struct Logger {
  quill::Logger *impl;
};

// Exposed because the CSV writers attach to the same sink.
std::shared_ptr<quill::Sink> get_console_sink();

} // namespace lbm::logging

// Quill's own macros already drop statements below
// QUILL_COMPILE_ACTIVE_LOG_LEVEL, so there is no compile-time check here.
// The null guard is what lets a logger that was never created be passed to a
// log statement without crashing.
#define LBM_DETAIL_QUILL_LOG(logger, quill_macro, ...)                         \
  do {                                                                         \
    if (::lbm::logging::Logger *lbm_logger_ = (logger)) {                      \
      quill_macro(lbm_logger_->impl, __VA_ARGS__);                             \
    }                                                                          \
  } while (false)

#define LBM_LOG_TRACE_L3(logger, ...)                                          \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_TRACE_L3, __VA_ARGS__)
#define LBM_LOG_TRACE_L2(logger, ...)                                          \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_TRACE_L2, __VA_ARGS__)
#define LBM_LOG_TRACE_L1(logger, ...)                                          \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_TRACE_L1, __VA_ARGS__)
#define LBM_LOG_DEBUG(logger, ...)                                             \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_DEBUG, __VA_ARGS__)
#define LBM_LOG_INFO(logger, ...)                                              \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_INFO, __VA_ARGS__)
#define LBM_LOG_NOTICE(logger, ...)                                            \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_NOTICE, __VA_ARGS__)
#define LBM_LOG_WARNING(logger, ...)                                           \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_WARNING, __VA_ARGS__)
#define LBM_LOG_ERROR(logger, ...)                                             \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_ERROR, __VA_ARGS__)
#define LBM_LOG_CRITICAL(logger, ...)                                          \
  LBM_DETAIL_QUILL_LOG(logger, QUILL_LOG_CRITICAL, __VA_ARGS__)
