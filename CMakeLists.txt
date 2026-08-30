# 3.21 is the real floor, not 3.18:
#   - target_sources() on an INTERFACE library needs 3.19;
#   - PROJECT_IS_TOP_LEVEL is only defined from 3.21 (with 3.18 the old
#     `if(${PROJECT_IS_TOP_LEVEL})` expanded to `if()` -> configure error).
cmake_minimum_required(VERSION 3.21)

project(lbm-simulation-lib LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Global toolchain settings belong to whoever drives the build. When this
# project is add_subdirectory()'d we must not overwrite the parent's choices;
# the standard is requested per-target below via target_compile_features().
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# ------------------------------ OPTIONS -------------------------------------

option(LBM_ENABLE_CUDA      "Build the CUDA backend"                        OFF)
option(LBM_ENABLE_CONFIG    "Build the TOML configuration support (toml++)" ON)
option(LBM_ENABLE_PROFILING "Enable chrono-based profiling instrumentation" OFF)

# LBM_ENABLE_LOGGING, LBM_LOG_BACKEND (quill | ostream | none) and
# LBM_LOG_LEVEL are declared in cmake/Logging.cmake, next to the code that
# reads them. LBM_SANITIZE is declared in cmake/Sanitizers.cmake.
option(LBM_BUILD_SIMULATIONS "Build the executables under simulations/"
       ${PROJECT_IS_TOP_LEVEL})
option(LBM_BUILD_DOCS       "Configure the docs/ target"
       ${PROJECT_IS_TOP_LEVEL})

# ----------------------------- SANITIZERS -----------------------------------
# Must come before include(Logging) and before the toml++ fetch below: the
# flags are added with add_compile_options() at directory scope, and a
# FetchContent'd subproject only inherits what was set before its
# add_subdirectory(). Quill is a static library that ends up in every
# executable, so it has to be instrumented identically or the link fails
# (MSVC) / the shadow bookkeeping disagrees (GCC, Clang).
#
# Reads LBM_ENABLE_CUDA, so it also has to come after the option() above.
include(Sanitizers)

# ------------------------------ LOGGING -------------------------------------
# Defines lbm::logging and resolves the backend. Only include/lbm/logging/quill/
# knows Quill exists; with any other backend nothing is fetched and no quill
# header is reachable from the build -- so the target below can link it
# unconditionally.
include(Logging)

# -------------------------------- CORE --------------------------------------

set(LBM_SIM_LIB lbm-sim)

file(GLOB_RECURSE LBM_CORE_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.hpp")

file(GLOB_RECURSE LBM_SHARED_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm/format/*.hpp")

add_library(${LBM_SIM_LIB} INTERFACE)
add_library(lbm::sim ALIAS ${LBM_SIM_LIB})

# PRIVATE on an INTERFACE library: the headers show up in IDEs/compile_commands
# but are not injected as sources into every consumer, which is what the old
# `add_library(... INTERFACE ${headers})` did.
target_sources(${LBM_SIM_LIB} PRIVATE ${LBM_CORE_HEADERS} ${LBM_SHARED_HEADERS})

target_include_directories(${LBM_SIM_LIB} INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_features(${LBM_SIM_LIB} INTERFACE cxx_std_17)

target_link_libraries(${LBM_SIM_LIB} INTERFACE lbm::logging)

# -------------------------------- OPENMP ------------------------------------

if(LBM_SANITIZE_DISABLES_OPENMP)
  # MSVC's AddressSanitizer does not support /openmp, so an ASan build on
  # Windows is necessarily serial. The `#pragma omp` directives stay in the
  # source and are simply ignored; silence the unknown-pragma warning so the
  # build output stays readable.
  message(WARNING
    "OpenMP is disabled because MSVC's AddressSanitizer does not support "
    "/openmp. The CPU solvers in this build run single-threaded -- fine for "
    "hunting memory bugs, useless for timings.")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:/wd4068>")
else()
  # MSVC's default OpenMP is 2.0, which has no unsigned loop counters and no
  # collapse. From CMake 3.30 the runtime is selectable here; before that,
  # FindOpenMP hands out plain /openmp and adding /openmp:experimental on top
  # of it is a hard cl.exe error (the two options are mutually exclusive),
  # which is what the old `target_compile_options(... /openmp:experimental)`
  # was doing.
  if(MSVC AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    set(OpenMP_RUNTIME_MSVC "experimental")
  endif()

  find_package(OpenMP REQUIRED)
  target_link_libraries(${LBM_SIM_LIB} INTERFACE OpenMP::OpenMP_CXX)
endif()

# ------------------------------- SETTINGS -----------------------------------

if(LBM_ENABLE_PROFILING)
  target_compile_definitions(${LBM_SIM_LIB} INTERFACE LBM_PROFILING)
endif()

# -Wsign-compare is a GCC/Clang C++ flag: it was reaching MSVC (which rejects
# it) and nvcc (which needs -Xcompiler for host flags). Fence it in.
target_compile_options(${LBM_SIM_LIB} INTERFACE
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-Wsign-compare>)

# ------------------------- CONFIG (part of the core) ------------------------
# TOML parsing is a feature of the library, not a separate target: toml++ and
# include/lbm/config/ are attached straight to lbm::sim. Consumers get it by
# linking lbm::sim -- there is nothing extra to link.
#
# With LBM_ENABLE_CONFIG=OFF nothing is fetched, include/lbm/config/ is not part
# of the build, LBM_ENABLE_CONFIG is not defined, and code that #includes those
# headers is expected not to be built either (`#if LBM_ENABLE_CONFIG`).
if(LBM_ENABLE_CONFIG)
  include(FetchContent)
  FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(tomlplusplus)

  target_link_libraries(${LBM_SIM_LIB} INTERFACE tomlplusplus::tomlplusplus)

  # For code that can do without a config file: #if LBM_ENABLE_CONFIG.
  target_compile_definitions(${LBM_SIM_LIB} INTERFACE LBM_ENABLE_CONFIG)
endif()

# --------------------------------- CUDA -------------------------------------
# The module no longer wraps itself in `if(LBM_ENABLE_CUDA)`: the guard lives
# here, so reading this file tells you what gets built. It must come before
# add_subdirectory(simulations) -- the CUDA language and
# CMAKE_CUDA_SEPARABLE_COMPILATION have to be in place before the cuda_*
# executables are created.
if(LBM_ENABLE_CUDA)
  include(Cuda)
endif()

# ------------------------------ SUBPROJECTS ---------------------------------

if(LBM_BUILD_SIMULATIONS)
  add_subdirectory(simulations)
endif()

if(LBM_BUILD_DOCS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/docs/CMakeLists.txt")
  add_subdirectory(docs)
endif()

# ------------------------------- SUMMARY ------------------------------------

message(STATUS "lbm-simulation-lib configuration:")
message(STATUS "  logging   : ${LBM_LOG_BACKEND} (${LBM_LOG_LEVEL})")
message(STATUS "  config    : ${LBM_ENABLE_CONFIG}")
message(STATUS "  cuda      : ${LBM_ENABLE_CUDA}")
message(STATUS "  profiling : ${LBM_ENABLE_PROFILING}")
if(LBM_SANITIZE_ACTIVE)
  message(STATUS "  sanitize  : ${LBM_SANITIZE}")
else()
  message(STATUS "  sanitize  : off")
endif()
