#include "lbm/format/format.hpp"
#include "lbm/logging/logging.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <unordered_map>

namespace lbm::logging {

namespace {

std::mutex &registry_mutex() {
  static std::mutex mutex;
  return mutex;
}

std::unordered_map<std::string, std::unique_ptr<Logger>> &registry() {
  static std::unordered_map<std::string, std::unique_ptr<Logger>> loggers;
  return loggers;
}

// Serialises whole lines: without it two threads interleave mid-message.
std::mutex &output_mutex() {
  static std::mutex mutex;
  return mutex;
}

std::tm local_time(std::time_t t) noexcept {
  std::tm tm{};
#ifdef _WIN32
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  return tm;
}

} // namespace

void setup() {}

void shutdown() {
  std::lock_guard<std::mutex> lock(output_mutex());
  std::clog.flush();
}

Logger *create_or_get_logger(std::string const &name) {
  std::lock_guard<std::mutex> lock(registry_mutex());

  if (auto it = registry().find(name); it != registry().end()) {
    return it->second.get();
  }

  auto handle = std::make_unique<Logger>(
      Logger{name, level_from_string(LBM_LOG_LEVEL_STR)});
  return registry().emplace(name, std::move(handle)).first->second.get();
}

void set_log_level(Logger *logger, LogLevel level) {
  if (logger != nullptr) {
    logger->level = level;
  }
}

namespace detail {

void write_line(Logger const &logger, char const *level_tag,
                std::string const &message) {
  auto const now = std::chrono::system_clock::now();
  auto const secs = std::chrono::system_clock::to_time_t(now);
  auto const ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      now.time_since_epoch()) %
                  std::chrono::seconds(1);
  std::tm const tm = local_time(secs);

  // La riga viene composta fuori dal lock: sotto il mutex resta solo la
  // scrittura, e un singolo operator<< non può essere spezzato a metà.
  std::ostringstream line;
  line << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << '.' << std::setfill('0')
       << std::setw(3) << ms.count() << " [" << logger.name << "] ["
       << level_tag << "] " << message << '\n';

  std::lock_guard<std::mutex> lock(output_mutex());
  std::clog << line.str();
}

} // namespace detail
} // namespace lbm::logging
