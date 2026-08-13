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

quill::Logger *create_or_get_logger(const std::string &logger_name) {
  return quill::Frontend::create_or_get_logger(logger_name, get_console_sink());
}

} // namespace logging
} // namespace lbm
