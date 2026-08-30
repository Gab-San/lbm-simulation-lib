/**
 * @file csv-writer.hpp
 * @brief CsvWriter: a synchronous, schema-driven CSV writer.
 */

#pragma once

/// Synchronous CSV writer, independent of LBM_LOG_BACKEND.
///
/// Takes the same schema shape as quill::CsvWriter - a type with static
/// `header` and `format` members - so existing schemas are reused unchanged:
///
///   struct ProfilingSchemaOpenMP {
///     static constexpr char const *header =
///         "size,collision_model,backend,n_threads,time";
///     static constexpr char const *format = "{},{},{},{},{:.2f}";
///   };
///
///   lbm::formatting::CsvWriter<ProfilingSchemaOpenMP> w{"profiling.csv"};
///   w.append_row(size, model, backend, n_threads, elapsed);
///
/// Rows are formatted on the calling thread and land in the file in call
/// order. Nothing is queued, so no row can be lost by skipping a shutdown or
/// by a backend thread that was never started, and the file does not depend on
/// how the project was configured.
///
/// Header-only, and part of lbm::sim rather than a target of its own: there is
/// nothing to compile and the include root is already public.

#include "lbm/format/format.hpp"

#include <filesystem>
#include <fstream>
#include <ios>
#include <string>

namespace lbm::format {

template <class Schema> class CsvWriter {
public:
  /// Opens `path` and writes the header row.
  ///
  /// With append=false the file is truncated. With append=true rows are added
  /// to whatever is there and the header is written only if the file was
  /// empty, so several processes -- a sweep that runs one executable per
  /// configuration -- can accumulate into one file.
  ///
  /// Throws std::ios_base::failure if the file cannot be opened: a run that
  /// cannot record its results has usually lost its point, and a silently
  /// missing results file is worse than a loud failure at startup.
  explicit CsvWriter(std::filesystem::path const &path, bool append = false)
      : out_(path, append ? (std::ios::out | std::ios::app) : std::ios::out) {
    if (!out_.is_open()) {
      throw std::ios_base::failure("lbm::format::CsvWriter: cannot open " +
                                   path.string());
    }
    // Opening created the file if it was missing, so this is the size of what
    // was already there. Nothing has been written through out_ yet.
    if (!append || std::filesystem::file_size(path) == 0) {
      out_ << Schema::header << '\n';
    }
  }

  CsvWriter(CsvWriter const &) = delete;
  CsvWriter &operator=(CsvWriter const &) = delete;

  CsvWriter(CsvWriter &&) noexcept = default;
  CsvWriter &operator=(CsvWriter &&) noexcept = default;

  ~CsvWriter() { close(); }

  /// Formats one row through Schema::format and terminates it with a newline.
  /// Arity is not checked: surplus arguments are ignored, missing ones expand
  /// to nothing. A no-op once close() has been called.
  template <class... Ts> void append_row(Ts const &...args) {
    if (!out_.is_open()) {
      return;
    }
    lbm::format::format_to(out_, Schema::format, args...);
    out_ << '\n';
  }

  void flush() {
    if (out_.is_open()) {
      out_.flush();
    }
  }

  /// Idempotent - safe to call explicitly and again from the destructor, which
  /// is what ProfilerWriter's destructor relies on.
  void close() {
    if (out_.is_open()) {
      out_.close();
    }
  }

  [[nodiscard]] bool is_open() const noexcept { return out_.is_open(); }

private:
  std::ofstream out_;
};

} // namespace lbm::format
