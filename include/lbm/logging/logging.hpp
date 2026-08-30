/**
 * @file logging.hpp
 * @brief Backend-agnostic logging facade.
 */

#pragma once

// Backend-agnostic logging facade.
//
// This is the only logging header the rest of the project includes. It never
// names a backend type: `Logger` is declared here and defined by whichever
// backend header is pulled in at the bottom of this file.
//
// Selected by CMake (cmake/Logging.cmake) through LBM_LOG_BACKEND.

#define LBM_LOG_BACKEND_QUILL 1
#define LBM_LOG_BACKEND_OSTREAM 2
#define LBM_LOG_BACKEND_NONE 3

#ifndef LBM_LOG_BACKEND
// Only reached when a header is compiled outside the build system.
#define LBM_LOG_BACKEND LBM_LOG_BACKEND_NONE
#endif

#include <string>
#include <string_view>

namespace lbm::logging {

// Ordered lowest to highest: backends filter with `logger_level <= statement`.
// The order matches LBM_LOG_LEVEL_VALUES in cmake/Logging.cmake, which is what
// makes LBM_ACTIVE_LOG_LEVEL below a plain list index.
enum class LogLevel {
  TraceL3 = 0,
  TraceL2 = 1,
  TraceL1 = 2,
  Debug = 3,
  Info = 4,
  Notice = 5,
  Warning = 6,
  Error = 7,
  Critical = 9,
  None = 10,
};

inline LogLevel level_from_string(std::string_view s) noexcept {
  if (s == "TRACE_L3") {
    return LogLevel::TraceL3;
  }
  if (s == "TRACE_L2") {
    return LogLevel::TraceL2;
  }
  if (s == "TRACE_L1") {
    return LogLevel::TraceL1;
  }
  if (s == "DEBUG") {
    return LogLevel::Debug;
  }
  if (s == "INFO") {
    return LogLevel::Info;
  }
  if (s == "NOTICE") {
    return LogLevel::Notice;
  }
  if (s == "WARNING") {
    return LogLevel::Warning;
  }
  if (s == "ERROR") {
    return LogLevel::Error;
  }
  if (s == "CRITICAL") {
    return LogLevel::Critical;
  }
  if (s == "NONE") {
    return LogLevel::None;
  }
  return LogLevel::Info; // fallback
}

#ifndef LBM_ACTIVE_LOG_LEVEL
#define LBM_ACTIVE_LOG_LEVEL 4 // Info
#endif

// Statements below this level are removed at compile time, before the runtime
// check against the logger's own level. Quill does this itself through
// QUILL_COMPILE_ACTIVE_LOG_LEVEL; the ostream backend uses this constant.
inline constexpr LogLevel compile_time_level =
    static_cast<LogLevel>(LBM_ACTIVE_LOG_LEVEL);

inline constexpr bool enabled_at_compile_time(LogLevel level) noexcept {
  return level >= compile_time_level;
}

/// Opaque handle. Core code holds a `Logger*` and passes it to the LBM_LOG_*
/// macros; it never touches the members, which differ per backend. Levels are
/// changed through set_log_level() below, never through a method.
struct Logger;

/// Starts the backend. Call once, at the top of main().
void setup();

/// Flushes and stops the backend. Call before leaving main() if the tail of
/// the log matters.
void shutdown();

/// Loggers are owned by the backend and live until shutdown(); the returned
/// pointer is stable. Calling twice with the same name returns the same
/// logger.
Logger *create_or_get_logger(std::string const &name);

/// Raises or lowers the run-time threshold of one logger. Statements below
/// the compile-time level are gone from the binary and cannot be re-enabled
/// here.
void set_log_level(Logger *logger, LogLevel level);

namespace detail {
// Never called at runtime: it exists so the discarding macros still
// type-check their arguments.
template <class... Ts> constexpr void ignore(Ts const &...) noexcept {}
} // namespace detail

} // namespace lbm::logging

// Each backend header defines `struct Logger` and the eight LBM_LOG_* macros.
#if LBM_LOG_BACKEND == LBM_LOG_BACKEND_QUILL
#include "lbm/logging/quill/backend.hpp"
#elif LBM_LOG_BACKEND == LBM_LOG_BACKEND_OSTREAM
#include "lbm/logging/ostream/backend.hpp"
#elif LBM_LOG_BACKEND == LBM_LOG_BACKEND_NONE
#include "lbm/logging/none/backend.hpp"
#else
#error "LBM_LOG_BACKEND must be LBM_LOG_BACKEND_{QUILL,OSTREAM,NONE}"
#endif
