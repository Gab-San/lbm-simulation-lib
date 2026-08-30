/**
 * @file backend.hpp
 * @brief The no-op logging backend.
 */

#pragma once

// No-op backend. Included by lbm/logging/logging.hpp; never include directly.
//
// `if (false)` and not a bare (void) cast: the arguments are still compiled,
// so a typo in a log statement is still an error and a variable used only in
// logging is still used. Nothing is evaluated at run time and the branch is
// folded away even at -O0.

namespace lbm::logging {
struct Logger {};
} // namespace lbm::logging

#define LBM_DETAIL_DISCARD_LOG(logger, ...)                                    \
  do {                                                                         \
    if (false) {                                                               \
      ::lbm::logging::detail::ignore(logger, __VA_ARGS__);                     \
    }                                                                          \
  } while (false)

#define LBM_LOG_TRACE_L3(logger, ...)                                          \
  LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_TRACE_L2(logger, ...)                                          \
  LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_TRACE_L1(logger, ...)                                          \
  LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_DEBUG(logger, ...) LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_INFO(logger, ...) LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_NOTICE(logger, ...) LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_WARNING(logger, ...) LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_ERROR(logger, ...) LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
#define LBM_LOG_CRITICAL(logger, ...)                                          \
  LBM_DETAIL_DISCARD_LOG(logger, __VA_ARGS__)
