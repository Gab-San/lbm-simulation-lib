/**
 * @file profiling.hpp
 * @brief Scope timers and the CSV they are dumped to.
 *
 * @c PROFILE_SCOPE(name) times the enclosing scope and accumulates total
 * time and call count per label in a process-wide registry; the solver dumps
 * the registry through a CsvWriter when the run ends.
 *
 * Active only when @c LBM_PROFILING is defined, which
 * @c -DLBM_ENABLE_PROFILING=ON does. Without it the macro expands to
 * nothing, so the instrumentation costs neither time nor code size in a
 * normal build.
 *
 * @see the "Output formats" page for the CSV schema.
 */

// lbm-sim/omp/profiling.hpp

#pragma once

/// One CSV writer per schema, shared by every call site in the process.
///
/// Owns the file's lifetime: open() once at startup, rows from anywhere,
/// close() at shutdown (or leave it to the destructor). Differences from the
/// obvious `static Writer` singleton, each of which is load-bearing:
///
///   - The writer is held in an optional, not a member of the singleton type.
///     CsvWriter has no default constructor - it opens a file or throws - so a
///     singleton that holds one by value cannot be constructed before a path
///     is known, which is exactly when the Meyers static runs.
///
///   - There is no getWriter(). Rows go through append_row(), so no reference
///     to the writer escapes: it cannot outlive a close(), and the mutex below
///     actually covers every write.
///
///   - append_row() is serialised. Rows are typically written once per run,
///     but an OpenMP region that records per-thread timings would otherwise
///     interleave two rows into one line, which corrupts the file rather than
///     just reordering it.
///
/// Cross-process sweeps: open(path, /*append=*/true) from each executable.
/// The header is written only into an empty file, so one CSV accumulates every
/// configuration. Within one process, append=false and one open() is simpler.

#include "lbm/format/csv-writer.hpp"

#include <filesystem>
#include <mutex>
#include <optional>
#include <utility>

namespace lbm {

struct ProfilingSchemaOpenMP {
  static constexpr char const *header =
      "label,size,collision_model,backend,n_threads,total,avg,calls";
  static constexpr char const *format = "{},{},{},{},{},{:.2f},{},{}";
};

namespace profiling {

template <class Schema> class Profiler {
public:
  static Profiler &get() {
    // Thread-safe initialisation since C++11; destroyed at exit, which flushes
    // and closes the file.
    static Profiler instance;
    return instance;
  }

  Profiler(Profiler const &) = delete;
  Profiler &operator=(Profiler const &) = delete;
  Profiler(Profiler &&) = delete;
  Profiler &operator=(Profiler &&) = delete;

  /// Opens the output file, replacing any file already open. Throws
  /// std::ios_base::failure if it cannot be opened.
  void open(std::filesystem::path const &path, bool append = false) {
    std::lock_guard<std::mutex> const lock(mutex_);
    writer_.reset(); // close the previous file before opening another
    std::filesystem::create_directories(path.parent_path());
    writer_.emplace(path, append);
  }

  /// Writes one row. A no-op if open() was never called - instrumentation left
  /// in a build that does not want a profile should cost nothing and say
  /// nothing, and the LBM_PROFILING macro is the switch that expresses intent.
  template <class... Ts> void append_row(Ts const &...args) {
    std::lock_guard<std::mutex> const lock(mutex_);
    if (writer_) {
      writer_->append_row(args...);
    }
  }

  void flush() {
    std::lock_guard<std::mutex> const lock(mutex_);
    if (writer_) {
      writer_->flush();
    }
  }

  /// Idempotent.
  void close() {
    std::lock_guard<std::mutex> const lock(mutex_);
    writer_.reset();
  }

  [[nodiscard]] bool is_open() const {
    std::lock_guard<std::mutex> const lock(mutex_);
    return writer_.has_value() && writer_->is_open();
  }

private:
  Profiler() = default;
  ~Profiler() = default;

  mutable std::mutex mutex_;
  std::optional<lbm::format::CsvWriter<Schema>> writer_;
};

} // namespace profiling
} // namespace lbm

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
