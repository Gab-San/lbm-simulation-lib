// lbm-sim/omp/profiling.hpp
#pragma once

#include <chrono>
#include <string>
#include <unordered_map>

#ifdef LBM_PROFILING

#include <fstream>
#include <iostream>
#include <mutex>

namespace lbm::profiling {

#define LBM_CONCAT_(a, b) a##b
#define LBM_CONCAT(a, b) LBM_CONCAT_(a, b)

struct TimingEntry {
  double total_ms = 0.0;
  std::size_t calls = 0;
};

inline std::unordered_map<std::string, TimingEntry> &registry() {
  static std::unordered_map<std::string, TimingEntry> reg;
  return reg;
}
inline std::mutex &registry_mutex() {
  static std::mutex m;
  return m;
}

class ScopedTimer {
public:
  explicit ScopedTimer(const char *label)
      : label_(label), start_(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    double ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - start_)
                    .count();
    std::lock_guard<std::mutex> lock(registry_mutex());
    auto &e = registry()[label_];
    e.total_ms += ms;
    e.calls += 1;
  }

private:
  const char *label_;
  std::chrono::steady_clock::time_point start_;
};

inline void reset() {
  std::lock_guard<std::mutex> lock(registry_mutex());
  registry().clear();
}

inline void report(std::ostream &os = std::cout) {
  std::lock_guard<std::mutex> lock(registry_mutex());
  for (auto const &[label, e] : registry()) {
    double avg = e.calls ? e.total_ms / static_cast<double>(e.calls) : 0.0;
    os << label << ": total=" << e.total_ms << "ms, calls=" << e.calls
       << ", avg=" << avg << "ms\n";
  }
}

inline void dump_csv(const std::string &path) {
  std::lock_guard<std::mutex> lock(registry_mutex());
  std::ofstream out(path);
  out << "label,total_ms,calls,avg_ms\n";
  for (auto const &[label, e] : registry()) {
    double avg = e.calls ? e.total_ms / static_cast<double>(e.calls) : 0.0;
    out << label << ',' << e.total_ms << ',' << e.calls << ',' << avg << '\n';
  }
}

} // namespace lbm::profiling

#define PROFILE_SCOPE(name)                                                    \
  ::lbm::profiling::ScopedTimer LBM_CONCAT(_prof_, __LINE__)(name)

#else

#define PROFILE_SCOPE(name)                                                    \
  do {                                                                         \
  } while (0)

#endif // LBM_PROFILING
