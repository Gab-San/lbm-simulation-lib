#include "lbm/logging/logging.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
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

bool is_digit(char c) noexcept { return c >= '0' && c <= '9'; }

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

namespace detail {

void apply_spec(std::ostream &os, std::string_view spec) {
  if (spec.empty()) {
    return;
  }

  if (auto const dot = spec.find('.'); dot != std::string_view::npos) {
    int precision = 0;
    for (std::size_t i = dot + 1; i < spec.size() && is_digit(spec[i]); ++i) {
      precision = precision * 10 + (spec[i] - '0');
    }
    os << std::setprecision(precision);
  }

  switch (spec.back()) {
  case 'f':
  case 'F':
    os << std::fixed;
    break;
  case 'e':
  case 'E':
    os << std::scientific;
    break;
  case 'g':
  case 'G':
    os.unsetf(std::ios_base::floatfield);
    break;
  case 'x':
    os << std::hex;
    break;
  case 'X':
    os << std::hex << std::uppercase;
    break;
  case 'o':
    os << std::oct;
    break;
  case 'd':
    os << std::dec;
    break;
  default:
    break; // width, fill and alignment are not supported
  }
}

void vformat_to(std::ostream &os, std::string_view fmt, Arg const *args,
                std::size_t count) {
  std::size_t next = 0;

  for (std::size_t i = 0; i < fmt.size(); ++i) {
    char const c = fmt[i];

    if (c == '{') {
      if (i + 1 < fmt.size() && fmt[i + 1] == '{') { // {{ -> {
        os << '{';
        ++i;
        continue;
      }

      auto const close = fmt.find('}', i + 1);
      if (close == std::string_view::npos) { // unterminated: emit as-is
        os << fmt.substr(i);
        return;
      }

      std::string_view spec = fmt.substr(i + 1, close - i - 1);
      if (!spec.empty() && spec.front() == ':') {
        spec.remove_prefix(1);
      }

      if (next < count) {
        args[next].print(os, args[next].value, spec);
        ++next;
      } else {
        os << "{?}"; // more placeholders than arguments
      }

      i = close;
      continue;
    }

    if (c == '}' && i + 1 < fmt.size() && fmt[i + 1] == '}') { // }} -> }
      os << '}';
      ++i;
      continue;
    }

    os << c;
  }
}

void write_line(Logger const &logger, char const *level_tag,
                std::string const &message) {
  auto const now = std::chrono::system_clock::now();
  auto const seconds = std::chrono::system_clock::to_time_t(now);
  auto const millis = std::chrono::duration_cast<std::chrono::milliseconds>(
                          now.time_since_epoch()) %
                      1000;
  std::tm const tm = local_time(seconds);

  // Built outside the lock; only the write is serialised.
  std::ostringstream line;
  line << std::put_time(&tm, "%H:%M:%S") << '.' << std::setfill('0')
       << std::setw(3) << millis.count() << " [" << level_tag << "] "
       << logger.name << " - " << message << '\n';

  std::lock_guard<std::mutex> lock(output_mutex());
  std::clog << line.str();
  std::clog.flush(); // a fallback backend that loses the last lines on a
                     // crash is worse than a slow one
}

} // namespace detail

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

} // namespace lbm::logging
