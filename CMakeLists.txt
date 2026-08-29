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
option(LBM_ENABLE_CONFIG    "Build the TOML configuration parser"           ON)
option(LBM_ENABLE_PROFILING "Enable chrono-based profiling instrumentation" OFF)

# LBM_ENABLE_LOGGING, LBM_LOG_BACKEND (quill | ostream | none) and
# LBM_LOG_LEVEL are declared in cmake/Logging.cmake, next to the code that
# reads them.
option(LBM_BUILD_SIMULATIONS "Build the executables under simulations/"
       ${PROJECT_IS_TOP_LEVEL})
option(LBM_BUILD_DOCS       "Configure the docs/ target"
       ${PROJECT_IS_TOP_LEVEL})

# ------------------------------ LOGGING -------------------------------------
# Defines lbm::logging and resolves the backend. Only include/lbm/logging/quill/
# knows Quill exists; with any other backend nothing is fetched and no quill
# header is reachable from the build -- so the target below can link it
# unconditionally.
include(Logging)

# ------------------------------- CONFIG -------------------------------------
# Defines lbm::config (toml++), or an empty stub when LBM_ENABLE_CONFIG=OFF.
# Deliberately NOT linked into lbm::sim -- see the note at the bottom of
# cmake/Config.cmake.
include(Config)

# -------------------------------- CORE --------------------------------------

set(LBM_SIM_LIB lbm-sim)

file(GLOB_RECURSE LBM_CORE_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.hpp")

# include/lbm-sim/config/ belongs to lbm::config, which lists it itself. Leaving
# it here too would mean the headers stay in the build (in IDEs, in the
# generated project) after LBM_ENABLE_CONFIG=OFF, which is exactly the coupling
# the split is meant to remove.
list(FILTER LBM_CORE_HEADERS EXCLUDE REGEX "/include/lbm-sim/config/")

add_library(${LBM_SIM_LIB} INTERFACE)
add_library(lbm::sim ALIAS ${LBM_SIM_LIB})

# PRIVATE on an INTERFACE library: the headers show up in IDEs/compile_commands
# but are not injected as sources into every consumer, which is what the old
# `add_library(... INTERFACE ${headers})` did.
target_sources(${LBM_SIM_LIB} PRIVATE ${LBM_CORE_HEADERS})

target_include_directories(${LBM_SIM_LIB} INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_features(${LBM_SIM_LIB} INTERFACE cxx_std_17)

# MSVC's default OpenMP is 2.0, which has no unsigned loop counters and no
# collapse. From CMake 3.30 the runtime is selectable here; before that,
# FindOpenMP hands out plain /openmp and adding /openmp:experimental on top of
# it is a hard cl.exe error (the two options are mutually exclusive), which is
# what the old `target_compile_options(... /openmp:experimental)` was doing.
if(MSVC AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
  set(OpenMP_RUNTIME_MSVC "experimental")
endif()

find_package(OpenMP REQUIRED)

target_link_libraries(${LBM_SIM_LIB} INTERFACE
  OpenMP::OpenMP_CXX
  lbm::logging)

if(LBM_ENABLE_PROFILING)
  target_compile_definitions(${LBM_SIM_LIB} INTERFACE LBM_PROFILING)
endif()

# -Wsign-compare is a GCC/Clang C++ flag: it was reaching MSVC (which rejects
# it) and nvcc (which needs -Xcompiler for host flags). Fence it in.
target_compile_options(${LBM_SIM_LIB} INTERFACE
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-Wsign-compare>)

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
