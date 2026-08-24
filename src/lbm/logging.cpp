#include "lbm/logging.hpp"

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/sinks/ConsoleSink.h"

namespace lbm {
namespace logging {
void setup_quill() {
  quill::Backend::start();
  auto console_sink = get_console_sink();
}

std::shared_ptr<quill::Sink> get_console_sink() {
  static auto console_sink =
      quill::Frontend::create_or_get_sink<quill::ConsoleSink>("console_sink");
  return console_sink;
}

quill::LogLevel level_from_string(const std::string &s) {
  if (s == "TRACE_L3")
    return quill::LogLevel::TraceL3;
  if (s == "TRACE_L2")
    return quill::LogLevel::TraceL2;
  if (s == "TRACE_L1")
    return quill::LogLevel::TraceL1;
  if (s == "DEBUG")
    return quill::LogLevel::Debug;
  if (s == "INFO")
    return quill::LogLevel::Info;
  if (s == "WARNING")
    return quill::LogLevel::Warning;
  if (s == "ERROR")
    return quill::LogLevel::Error;
  if (s == "CRITICAL")
    return quill::LogLevel::Critical;
  return quill::LogLevel::Info; // fallback
}

quill::Logger *create_or_get_logger(const std::string &logger_name) {
  quill::Logger *logger =
      quill::Frontend::create_or_get_logger(logger_name, get_console_sink());
  logger->set_log_level(level_from_string(LBM_LOG_LEVEL_STR));
  return logger;
}

} // namespace logging
} // namespace lbm
