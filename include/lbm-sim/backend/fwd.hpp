/**
 * @file fwd.hpp
 * @brief Declaration-only view of BackendProperties.
 *
 * Include this to *name* the template; include backend/properties.hpp to
 * configure a backend. The split keeps @c <omp.h> and the CUDA runtime out
 * of headers that only mention the type in a signature.
 */

#pragma once

#include "lbm-sim/metadata.hpp"

namespace lbm::profiling {

/// Process-wide handle on the tuning knobs of an execution backend.
///
/// Only declared here: every backend provides its own specialization, so that
/// asking for one that has none fails at compile time instead of silently
/// yielding an empty object.
template <enum ExecutionBackend backend_t> class BackendProperties;

} // namespace lbm::profiling
