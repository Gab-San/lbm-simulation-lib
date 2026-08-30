#include "lbm/logging/logging.hpp"

// The LBM_LOG_* macros already expand to nothing in this configuration, so all
// that is left is keeping the four entry points linkable.

namespace lbm::logging {

void setup() {}

void shutdown() {}

Logger *create_or_get_logger(std::string const &name) {
  (void)name;
  static Logger null_logger{};
  return &null_logger;
}

void set_log_level(Logger *logger, LogLevel level) {
  (void)logger;
  (void)level;
}

} // namespace lbm::logging
