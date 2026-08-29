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

// --- fix per il bug di token-pasting con __LINE__ ---
// __LINE__ non viene espanso prima di ## senza un livello di indirezione:
// senza questa coppia di macro, PROFILE_SCOPE genera sempre la variabile
// letterale "_prof___LINE__" (redefinition se usata due volte nello stesso scope).
#define LBM_CONCAT_(a, b) a##b
#define LBM_CONCAT(a, b) LBM_CONCAT_(a, b)

struct TimingEntry {
  double total_ms = 0.0;
  std::size_t calls = 0;
};

inline std::unordered_map<std::string, TimingEntry>& registry() {
  static std::unordered_map<std::string, TimingEntry> reg;
  return reg;
}

// Protegge registry() da accessi concorrenti se PROFILE_SCOPE viene usato
// dentro una regione #pragma omp parallel (es. per profilare il kernel di
// collision per-thread). A livello di granularità "per fase" (collision,
// streaming, BC, halo) la contesa sul mutex è trascurabile; se in futuro
// si profila a grana più fine (es. per nodo dentro il loop) conviene
// passare a registry thread_local con merge esplicito a fine run.
inline std::mutex& registry_mutex() {
  static std::mutex m;
  return m;
}

class ScopedTimer {
public:
  explicit ScopedTimer(const char* label) : label_(label),
      start_(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start_).count();
    std::lock_guard<std::mutex> lock(registry_mutex());
    auto& e = registry()[label_];
    e.total_ms += ms;
    e.calls += 1;
  }
private:
  const char* label_;
  std::chrono::steady_clock::time_point start_;
};

inline void reset() {
  std::lock_guard<std::mutex> lock(registry_mutex());
  registry().clear();
}

inline void report(std::ostream& os = std::cout) {
  std::lock_guard<std::mutex> lock(registry_mutex());
  for (auto const& [label, e] : registry()) {
    double avg = e.calls ? e.total_ms / static_cast<double>(e.calls) : 0.0;
    os << label << ": total=" << e.total_ms << "ms, calls=" << e.calls
       << ", avg=" << avg << "ms\n";
  }
}

inline void dump_csv(const std::string& path) {
  std::lock_guard<std::mutex> lock(registry_mutex());
  std::ofstream out(path);
  out << "label,total_ms,calls,avg_ms\n";
  for (auto const& [label, e] : registry()) {
    double avg = e.calls ? e.total_ms / static_cast<double>(e.calls) : 0.0;
    out << label << ',' << e.total_ms << ',' << e.calls << ',' << avg << '\n';
  }
}

/* ------------------------------------------------------------------------
 * SAMPLE per estensione futura MPI (non funzionante, solo da adattare):
 * riduce le entry di registry() su rank 0 con min/max/avg per label,
 * per capire il load imbalance tra rank invece del solo tempo locale.
 *
 * inline void report_mpi(int rank, int nranks, std::ostream& os = std::cout) {
 *   for (auto const& [label, e] : registry()) {
 *     double local = e.total_ms;
 *     double sum, min, max;
 *     MPI_Reduce(&local, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 *     MPI_Reduce(&local, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
 *     MPI_Reduce(&local, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 *     if (rank == 0) {
 *       os << label << ": avg=" << (sum / nranks)
 *          << "ms, min=" << min << "ms, max=" << max
 *          << "ms (imbalance=" << (max - min) << "ms)\n";
 *     }
 *   }
 * }
 * ------------------------------------------------------------------------ */

} // namespace lbm::profiling

#define PROFILE_SCOPE(name) \
  ::lbm::profiling::ScopedTimer LBM_CONCAT(_prof_, __LINE__)(name)

#else

#define PROFILE_SCOPE(name) do {} while (0)

#endif // LBM_PROFILING