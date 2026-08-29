# lbm::logging -- the project's logging facade.
#
# The target always exists and is always linked into lbm::sim. What changes is
# which backend is compiled into it:
#
#   LBM_LOG_BACKEND=quill    include/lbm/logging/quill/    (fetches Quill)
#   LBM_LOG_BACKEND=ostream  include/lbm/logging/ostream/  (no dependency)
#   LBM_LOG_BACKEND=none     include/lbm/logging/none/     (macros expand away)
#
# Only include/lbm/logging/quill/ mentions quill; with any other backend
# nothing is fetched and no quill header is reachable from the build. That is
# why linking lbm::logging unconditionally is safe -- it is the *backend*, not
# the target, that carries the dependency.

include_guard(GLOBAL)

set(LBM_LOGGING_LIB lbm-logging)

# ------------------------------ SELECTION -----------------------------------

set(LBM_LOG_BACKEND "quill" CACHE STRING
    "Logging backend: quill | ostream | none")
set_property(CACHE LBM_LOG_BACKEND PROPERTY STRINGS quill ostream none)

# Convenience switch kept from the previous layout: OFF wins over the backend.
# A normal variable, so the cached choice is remembered if it is turned back on.
option(LBM_ENABLE_LOGGING "Emit log output (OFF forces LBM_LOG_BACKEND=none)" ON)
if(NOT LBM_ENABLE_LOGGING)
  set(LBM_LOG_BACKEND "none")
endif()

if(NOT LBM_LOG_BACKEND MATCHES "^(quill|ostream|none)$")
  message(FATAL_ERROR
    "LBM_LOG_BACKEND must be quill, ostream or none (got: ${LBM_LOG_BACKEND})")
endif()

# --------------------------- COMPILE-TIME LEVEL ------------------------------
#
# The order of this list is the numeric value of lbm::logging::LogLevel, which
# is what makes the list index below a valid LBM_ACTIVE_LOG_LEVEL. Keep the two
# in sync.

set(LBM_LOG_LEVEL_VALUES
    TRACE_L3 TRACE_L2 TRACE_L1 DEBUG INFO WARNING ERROR CRITICAL)

set(LBM_LOG_LEVEL "INFO" CACHE STRING "Compile-time active log level")
set_property(CACHE LBM_LOG_LEVEL PROPERTY STRINGS ${LBM_LOG_LEVEL_VALUES})

if(NOT LBM_LOG_LEVEL IN_LIST LBM_LOG_LEVEL_VALUES)
  message(FATAL_ERROR
    "LBM_LOG_LEVEL must be one of: ${LBM_LOG_LEVEL_VALUES} "
    "(got: ${LBM_LOG_LEVEL})")
endif()

list(FIND LBM_LOG_LEVEL_VALUES ${LBM_LOG_LEVEL} LBM_ACTIVE_LOG_LEVEL)

# ------------------------------- SOURCES -------------------------------------
#
# One directory per backend, in both trees, so the paths are derived from the
# backend name instead of an if/elseif chain. The only thing left to branch on
# is quill's extra header and its dependency.

set(LBM_LOGGING_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include/lbm/logging")
set(LBM_LOGGING_SRC_DIR     "${PROJECT_SOURCE_DIR}/src/lbm/logging")

set(LBM_LOGGING_SOURCES
    "${LBM_LOGGING_SRC_DIR}/${LBM_LOG_BACKEND}/backend.cpp")

set(LBM_LOGGING_HEADERS
    "${LBM_LOGGING_INCLUDE_DIR}/logging.hpp"
    "${LBM_LOGGING_INCLUDE_DIR}/${LBM_LOG_BACKEND}/backend.hpp")

if(LBM_LOG_BACKEND STREQUAL "quill")
  list(APPEND LBM_LOGGING_HEADERS
       "${LBM_LOGGING_INCLUDE_DIR}/quill/nvcc-compat.hpp")
endif()

# Cheap, because a missing backend .cpp otherwise surfaces as a pile of
# undefined references at link time in whichever executable happens to be first.
foreach(LBM_LOGGING_FILE IN LISTS LBM_LOGGING_SOURCES LBM_LOGGING_HEADERS)
  if(NOT EXISTS "${LBM_LOGGING_FILE}")
    message(FATAL_ERROR
      "LBM_LOG_BACKEND=${LBM_LOG_BACKEND} but this file is missing:\n"
      "  ${LBM_LOGGING_FILE}")
  endif()
endforeach()

# ----------------------------- DEPENDENCIES ----------------------------------

if(LBM_LOG_BACKEND STREQUAL "quill")
  include(FetchContent)

  # Set before MakeAvailable so quill's own option() calls pick them up; plain
  # CACHE without FORCE, which works whether or not quill's policy scope has
  # CMP0077 set to NEW.
  set(QUILL_BUILD_EXAMPLES   OFF CACHE BOOL "")
  set(QUILL_BUILD_TESTS      OFF CACHE BOOL "")
  set(QUILL_BUILD_BENCHMARKS OFF CACHE BOOL "")

  FetchContent_Declare(
    quill
    GIT_REPOSITORY https://github.com/odygrd/quill
    GIT_TAG        v12.1.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(quill)
endif()

# ------------------------------- TARGET --------------------------------------

add_library(${LBM_LOGGING_LIB} STATIC
            ${LBM_LOGGING_SOURCES} ${LBM_LOGGING_HEADERS})
add_library(lbm::logging ALIAS ${LBM_LOGGING_LIB})

target_include_directories(${LBM_LOGGING_LIB} PUBLIC
                           "${PROJECT_SOURCE_DIR}/include")

target_compile_features(${LBM_LOGGING_LIB} PUBLIC cxx_std_17)

# The rest of the project is header-only, so this is the one archive everything
# links. Keep it usable from a shared consumer.
set_target_properties(${LBM_LOGGING_LIB} PROPERTIES
                      POSITION_INDEPENDENT_CODE ON)

string(TOUPPER ${LBM_LOG_BACKEND} LBM_LOG_BACKEND_UPPER)
target_compile_definitions(${LBM_LOGGING_LIB} PUBLIC
  LBM_LOG_BACKEND=LBM_LOG_BACKEND_${LBM_LOG_BACKEND_UPPER}
  LBM_ACTIVE_LOG_LEVEL=${LBM_ACTIVE_LOG_LEVEL}
  LBM_LOG_LEVEL_STR="${LBM_LOG_LEVEL}")

if(LBM_LOG_BACKEND STREQUAL "quill")
  # PUBLIC and not PRIVATE on purpose: LBM_LOG_* expand to QUILL_LOG_* in the
  # consumer's translation unit, so consumers do need the quill headers. That
  # is the one place the dependency is genuinely public -- and it exists only
  # in this branch.
  target_link_libraries(${LBM_LOGGING_LIB} PUBLIC quill::quill)

  target_compile_definitions(${LBM_LOGGING_LIB} PUBLIC
    QUILL_COMPILE_ACTIVE_LOG_LEVEL=QUILL_COMPILE_ACTIVE_LOG_LEVEL_${LBM_LOG_LEVEL})

  # Quill's headers trip a couple of cudafe diagnostics when pulled into a .cu.
  # The suppression belongs to whoever brings in those headers, so it
  # disappears together with quill. SHELL: keeps "-Xcudafe <arg>" as two shell
  # words; without it CMake quotes the whole string as one argument.
  target_compile_options(${LBM_LOGGING_LIB} INTERFACE
    "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcudafe --diag_suppress=128,--diag_suppress=2417>")
endif()

message(STATUS "LBM logging: backend=${LBM_LOG_BACKEND} level=${LBM_LOG_LEVEL}")

# The profiling writer is still a quill::CsvWriter, so it needs that backend.
# Rejected loudly rather than silently dropping a feature that was asked for --
# and the honest fix is the ofstream CsvWriter fallback, after which this block
# goes away.
if(LBM_ENABLE_PROFILING AND NOT LBM_LOG_BACKEND STREQUAL "quill")
  message(FATAL_ERROR
    "LBM_ENABLE_PROFILING=ON currently requires LBM_LOG_BACKEND=quill "
    "(got: ${LBM_LOG_BACKEND}): lbm::profiling writes through quill::CsvWriter.")
endif()
