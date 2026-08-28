#pragma once

// FIXME: REWIRE
#include "lbm-sim/backend.hpp"

#include "lbm-sim/metadata.hpp"

namespace lbm::profiling {

/// Process-wide handle on the tuning knobs of an execution backend.
///
/// Only declared here: every backend provides its own specialization, so that
/// asking for one that has none fails at compile time instead of silently
/// yielding an empty object.
template <enum ExecutionBackend backend_t> class BackendProperties;

} // namespace lbm::profiling
