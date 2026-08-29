#include "lbm/logging/logging.hpp"

#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/sinks/ConsoleSink.h"

#include <memory>
#include <mutex>
#include <unordered_map>

namespace lbm::logging {

namespace {

quill::LogLevel to_quill(LogLevel level) noexcept {
  switch (level) {
  case LogLevel::TraceL3:
    return quill::LogLevel::TraceL3;
  case LogLevel::TraceL2:
    return quill::LogLevel::TraceL2;
  case LogLevel::TraceL1:
    return quill::LogLevel::TraceL1;
  case LogLevel::Debug:
    return quill::LogLevel::Debug;
  case LogLevel::Info:
    return quill::LogLevel::Info;
  case LogLevel::Notice:
    return quill::LogLevel::Notice;
  case LogLevel::Warning:
    return quill::LogLevel::Warning;
  case LogLevel::Error:
    return quill::LogLevel::Error;
  case LogLevel::Critical:
    return quill::LogLevel::Critical;
  case LogLevel::None:
    return quill::LogLevel::None;
  }
  return quill::LogLevel::Info;
}

// The handles are ours, so we own them. Node-based on purpose: the pointers
// handed out must stay valid as the map grows. Both are function-local so
// they are initialised on first use, whatever the static init order is.
std::mutex &registry_mutex() {
  static std::mutex mutex;
  return mutex;
}

std::unordered_map<std::string, std::unique_ptr<Logger>> &registry() {
  static std::unordered_map<std::string, std::unique_ptr<Logger>> loggers;
  return loggers;
}

} // namespace

std::shared_ptr<quill::Sink> get_console_sink() {
  static auto console_sink =
      quill::Frontend::create_or_get_sink<quill::ConsoleSink>("console_sink");
  return console_sink;
}

void setup() {
  quill::Backend::start();
  (void)get_console_sink();
}

void shutdown() { quill::Backend::stop(); }

Logger *create_or_get_logger(std::string const &name) {
  std::lock_guard<std::mutex> lock(registry_mutex());

  if (auto it = registry().find(name); it != registry().end()) {
    return it->second.get();
  }

  quill::Logger *impl =
      quill::Frontend::create_or_get_logger(name, get_console_sink());
  impl->set_log_level(to_quill(level_from_string(LBM_LOG_LEVEL_STR)));

  auto handle = std::make_unique<Logger>(Logger{impl});
  return registry().emplace(name, std::move(handle)).first->second.get();
}

void set_log_level(Logger *logger, LogLevel level) {
  if (logger != nullptr && logger->impl != nullptr) {
    logger->impl->set_log_level(to_quill(level));
  }
}

} // namespace lbm::logging
