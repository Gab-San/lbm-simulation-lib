// #pragma once

// #include "lbm-sim/backend/fwd.hpp"
// #include <omp.h>

// namespace lbm::profiling {

// /// This class stores the values to use, and allows them to be applied to the runtime.
// #if defined(_MSC_VER) && !defined(_OPENMP_LLVM_RUNTIME)
// using omp_schedule_type = int;
// #define LBM_HAS_OMP_SCHEDULE_API 0
// #else
// using omp_schedule_type = omp_sched_t;
// #define LBM_HAS_OMP_SCHEDULE_API 1
// #endif

// /// OpenMP knobs. Meant to be filled once at start-up from the configuration
// /// and only read afterwards: the accessors do no locking, so mutating them
// /// while a solver runs is a data race.
// template <> class BackendProperties<OPEN_MP> {
// public:
//   /// The single instance. The static local is initialized on first use, and
//   /// since C++11 that initialization is thread safe; being a local of an
//   /// implicitly inline function, it is also the same object in every
//   /// translation unit.
//   static BackendProperties &get() {
//     static BackendProperties instance;
//     return instance;
//   }

//   // A singleton must not be duplicated: a copy would be an object whose
//   // mutations nobody else can observe.
//   BackendProperties(const BackendProperties &) = delete;
//   BackendProperties &operator=(const BackendProperties &) = delete;
//   BackendProperties(BackendProperties &&) = delete;
//   BackendProperties &operator=(BackendProperties &&) = delete;

//   unsigned int getNumThreads() const noexcept { return num_threads; }
//   omp_schedule_type getSchedule() const noexcept { return schedule; }
//   int getChunkSize() const noexcept { return chunk_size; }
//   bool getDynamicThreads() const noexcept { return dynamic_threads; }

//   void setNumThreads(const unsigned int num_threads_) noexcept {
//     num_threads = num_threads_;
//   }

//   /// `chunk_size_` <= 0 leaves the chunk size to the runtime.
//   void setSchedule(const omp_schedule_type schedule_,
//                    const int chunk_size_ = 0) noexcept {
//     schedule = schedule_;
//     chunk_size = chunk_size_;
//   }

// //   /// With dynamic adjustment on, the runtime is free to hand a parallel region
// //   /// *fewer* threads than asked for, which makes a timing run unreproducible.
// //   /// Defaults to off; turn it back on only to measure that behaviour on
// //   /// purpose.
// //   void setDynamicThreads(const bool dynamic_threads_) noexcept {
// //     dynamic_threads = dynamic_threads_;
// //   }

//   /// Push the stored values into the OpenMP runtime.
//   ///
//   /// These are ICVs of the calling thread, so call this from the initial
//   /// thread and *outside* any parallel region: called inside one, it would
//   /// only affect the regions that thread goes on to nest. The schedule is
//   /// only honoured by loops declared `schedule(runtime)`.
//   void apply() const noexcept {
//     omp_set_dynamic(dynamic_threads ? 1 : 0);
//     omp_set_num_threads(static_cast<int>(num_threads));
// #if LBM_HAS_OMP_SCHEDULE_API
//     omp_set_schedule(schedule, chunk_size);
// #endif
//   }

// //   /// RAII counterpart of `apply()`: snapshots the runtime's ICVs, applies the
// //   /// stored properties, and restores the snapshot on scope exit. Use it to
// //   /// keep a profiled region from leaking its settings into whatever runs
// //   /// next:
// //   ///
// //   ///     {
// //   ///       const auto scope = BackendProperties<OPEN_MP>::get().scopedApply();
// //   ///       simulation.solve(solver);
// //   ///     } // runtime is back to what it was
// //   ///
// //   /// Same threading rule as `apply()`: construct it on the initial thread,
// //   /// outside any parallel region.
// //   class [[nodiscard]] Scope {
// //   public:
// //     explicit Scope(const BackendProperties &props) noexcept
// //         // omp_get_max_threads() reports the *nthreads-var* ICV, i.e. exactly
// //         // what omp_set_num_threads() last set, not a thread count in flight.
// //         : prev_num_threads(omp_get_max_threads()),
// //           prev_dynamic(omp_get_dynamic()) {
// //       omp_get_schedule(&prev_schedule, &prev_chunk_size);
// //       props.apply();
// //     }

//     ~Scope() {
//     #if LBM_HAS_OMP_SCHEDULE_API
//       omp_set_schedule(prev_schedule, prev_chunk_size);
//     #endif
//       omp_set_num_threads(prev_num_threads);
//       omp_set_dynamic(prev_dynamic);
//     }

//     // Restoring twice, or from the wrong scope, would defeat the purpose.
//     Scope(const Scope &) = delete;
//     Scope &operator=(const Scope &) = delete;
//     Scope(Scope &&) = delete;
//     Scope &operator=(Scope &&) = delete;

//   private:
//     int prev_num_threads;
//     int prev_dynamic;
//     omp_schedule_type prev_schedule{};
//     int prev_chunk_size;
//   };

//   /// Returned by value: C++17 guaranteed copy elision builds the `Scope`
//   /// straight into the caller's variable, so the deleted move is never
//   /// needed.
//   [[nodiscard]] Scope scopedApply() const { return Scope(*this); }

// private:
//   /// Defaults mirror what the runtime would have done on its own, so an
//   /// unconfigured process behaves exactly as before -- except for
//   /// `dynamic_threads`, which starts off because reproducible timing is the
//   /// whole point of applying these at all.
//   BackendProperties() noexcept
//       : num_threads(static_cast<unsigned int>(omp_get_max_threads())),
//         schedule(0), chunk_size(0), dynamic_threads(false) {}

//   ~BackendProperties() = default;

//   unsigned int num_threads;
//   omp_schedule_type schedule;
//   int chunk_size;
//   bool dynamic_threads;
// };

// } // namespace lbm::profiling

// #undef LBM_HAS_OMP_SCHEDULE_API
