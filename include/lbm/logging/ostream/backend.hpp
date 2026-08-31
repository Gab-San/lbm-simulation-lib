#pragma once

// std::ostream backend. Included by lbm/logging/logging.hpp; never include
// directly.
//
// Synchronous, formats on the calling thread, writes to std::clog under a
// mutex. It is the dependency-free fallback, not a replacement for quill on a
// hot path.
//
// Formatting is a small brace-substitution pass over operator<<, so that call
// sites written against quill keep working unchanged. What is supported:
//
//   {}          the argument, streamed with operator<<
//   {{  }}      literal braces
//   {:.3f}      precision + fixed / scientific / general  (f e g F E G)
//   {:x} {:X} {:o} {:d}   integer base
//
// Anything else in a spec is ignored and the argument is streamed plainly.
// A type is loggable here as long as it has an operator<<, which is what
// lbm::utils::Point and Vector already rely on through ostream_formatter.

#include <cstddef>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>

namespace lbm::logging {

struct Logger {
  std::string name;
  LogLevel level;
};

namespace detail {

// Type-erased argument: a pointer to the caller's value plus the function
// that knows how to stream it. No allocation, no std::function.
using ArgPrinter = void (*)(std::ostream &, void const *, std::string_view);

struct Arg {
  void const *value;
  ArgPrinter print;
};

// Applies what of the format spec we understand to the stream.
void apply_spec(std::ostream &os, std::string_view spec);

template <class T>
void print_arg(std::ostream &os, void const *value, std::string_view spec) {
  auto const flags = os.flags();
  auto const precision = os.precision();
  apply_spec(os, spec);
  os << *static_cast<T const *>(value);
  os.precision(precision);
  os.flags(flags);
}

template <class T> Arg make_arg(T const &value) noexcept {
  return Arg{static_cast<void const *>(std::addressof(value)), &print_arg<T>};
}

void vformat_to(std::ostream &os, std::string_view fmt, Arg const *args,
                std::size_t count);

void write_line(Logger const &logger, char const *level_tag,
                std::string const &message);

template <class... Ts>
void log_line(Logger const &logger, char const *level_tag, std::string_view fmt,
              Ts const &...args) {
  // Trailing sentinel: a zero-sized array is ill-formed when there are no
  // arguments beyond the format string.
  Arg const packed[] = {make_arg(args)..., Arg{nullptr, nullptr}};
  std::ostringstream out;
  vformat_to(out, fmt, packed, sizeof...(Ts));
  write_line(logger, level_tag, out.str());
}

} // namespace detail
} // namespace lbm::logging

// The parameters carry the lbm_ prefix on purpose: a macro parameter named
// `level` would also be substituted into `lbm_logger_->level`.
#define LBM_DETAIL_OSTREAM_LOG(lbm_logger_arg, lbm_level_arg, lbm_tag_arg,     \
                               ...)                                            \
  do {                                                                         \
    if constexpr (::lbm::logging::enabled_at_compile_time(lbm_level_arg)) {    \
      ::lbm::logging::Logger *lbm_logger_ = (lbm_logger_arg);                  \
      if (lbm_logger_ != nullptr && lbm_logger_->level <= (lbm_level_arg)) {   \
        ::lbm::logging::detail::log_line(*lbm_logger_, lbm_tag_arg,            \
                                         __VA_ARGS__);                         \
      }                                                                        \
    }                                                                          \
  } while (false)

#define LBM_LOG_TRACE_L3(logger, ...)                                          \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::TraceL3,            \
                         "TRACE_L3", __VA_ARGS__)
#define LBM_LOG_TRACE_L2(logger, ...)                                          \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::TraceL2,            \
                         "TRACE_L2", __VA_ARGS__)
#define LBM_LOG_TRACE_L1(logger, ...)                                          \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::TraceL1,            \
                         "TRACE_L1", __VA_ARGS__)
#define LBM_LOG_DEBUG(logger, ...)                                             \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Debug, "DEBUG",     \
                         __VA_ARGS__)
#define LBM_LOG_INFO(logger, ...)                                              \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Info, "INFO",       \
                         __VA_ARGS__)
#define LBM_LOG_NOTICE(logger, ...)                                            \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Notice, "NOTICE",   \
                         __VA_ARGS__)
#define LBM_LOG_WARNING(logger, ...)                                           \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Warning, "WARNING", \
                         __VA_ARGS__)
#define LBM_LOG_ERROR(logger, ...)                                             \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Error, "ERROR",     \
                         __VA_ARGS__)
#define LBM_LOG_CRITICAL(logger, ...)                                          \
  LBM_DETAIL_OSTREAM_LOG(logger, ::lbm::logging::LogLevel::Critical,           \
                         "CRITICAL", __VA_ARGS__)
