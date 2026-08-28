// lbm-sim/omp/annotations.hpp (aggiunta) o nuovo lbm-sim/profiling.hpp
#include <chrono>
#include <string>
#include <unordered_map>

#ifdef LBM_PROFILING

namespace lbm::profiling {

struct TimingEntry {
  double total_ms = 0.0;
  std::size_t calls = 0;
};

inline std::unordered_map<std::string, TimingEntry>& registry() {
  static std::unordered_map<std::string, TimingEntry> reg;
  return reg;
}

class ScopedTimer {
public:
  explicit ScopedTimer(const char* label) : label_(label),
      start_(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start_).count();
    auto& e = registry()[label_];
    e.total_ms += ms;
    e.calls += 1;
  }
private:
  const char* label_;
  std::chrono::steady_clock::time_point start_;
};

} // namespace lbm::profiling

#define PROFILE_SCOPE(name) ::lbm::profiling::ScopedTimer _prof_##__LINE__(name)

#else
#define PROFILE_SCOPE(name) do {} while (0)
#endif