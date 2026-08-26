cmake_minimum_required(VERSION 3.18)

project(lbm-simulation-lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(LBM_ENABLE_CUDA "Enable and build with CUDA Backend" OFF)

# ------------ QUILL LOGGING LIBRARY ------------
include(Quill)

set(LBM_LOGGING_HELPER_LIB lbm-logging)

# Solo i sorgenti del logging: una GLOB_RECURSE su src/lbm/*.cpp si tirava
# dentro anche src/lbm/config/config_parser.cpp, il vecchio parser del
# formato di configurazione custom, che e' orfano (nessuno lo usa) e non
# compila piu' contro l'API attuale. E' sostituito dal parser TOML in
# include/lbm-sim/config/.
set(
  LBM_LOGGING_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm/logging.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/lbm/logging.cpp"
)

add_library(${LBM_LOGGING_HELPER_LIB} STATIC ${LBM_LOGGING_FILES})
target_include_directories(
  ${LBM_LOGGING_HELPER_LIB}
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(${LBM_LOGGING_HELPER_LIB} PUBLIC quill::quill)
lbm_apply_quill_log_level(${LBM_LOGGING_HELPER_LIB})

# -------------------------------------------

# ------------ TOML PARSING LIBRARY ------------
include(Toml)
# -------------------------------------------

set(LBM_SIM_LIB lbm-sim)

file(
  GLOB_RECURSE
  LBM_CORE_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.hpp"
)

add_library(${LBM_SIM_LIB} INTERFACE ${LBM_CORE_HEADERS})

target_include_directories(
  ${LBM_SIM_LIB}
  INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

find_package(OpenMP REQUIRED)

target_link_libraries(${LBM_SIM_LIB} INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(${LBM_SIM_LIB} INTERFACE ${LBM_LOGGING_HELPER_LIB})
target_link_libraries(${LBM_SIM_LIB} INTERFACE tomlplusplus::tomlplusplus)

target_compile_options(${LBM_SIM_LIB} INTERFACE
    $<$<CONFIG:Debug>:-Wsign-compare>
)
# ------------------------------------------------------

if(MSVC)
  target_compile_options(${LBM_SIM_LIB} INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/openmp:experimental>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:/O2>
    )
endif()

# ------------------------------------------------------

include(Cuda)

# ------------------------------------------------------

if (${PROJECT_IS_TOP_LEVEL})
  add_subdirectory(simulations)
endif()

# ------------------------------------------------------

add_subdirectory(docs)
