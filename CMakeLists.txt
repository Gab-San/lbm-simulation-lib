cmake_minimum_required(VERSION 3.18)

project(lbm-simulation-lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

option(LBM_ENABLE_CUDA "Enable and build with CUDA Backend" OFF)

if (LBM_ENABLE_CUDA)
  enable_language(CUDA)
endif()

find_package(OpenMP REQUIRED)

file(
  GLOB_RECURSE
  LBM_CORE_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.hpp"
)

add_library(lib-lbm-sim INTERFACE ${LBM_CORE_HEADERS})
target_include_directories(
  lib-lbm-sim
  INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(lib-lbm-sim INTERFACE OpenMP::OpenMP_CXX)

if(MSVC)
    target_compile_options(lib-lbm-sim INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/openmp:experimental>
        # flag di ottimizzazione per release, sostituisce
        # il loop unrolling di OpenMP che non funziona con MSVC!!!
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:/O2>
    )
endif()

if(LBM_ENABLE_CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED True)

  # TODO: understand this!
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)   # adjust to your target GPU(s)
  endif()

  file(GLOB_RECURSE LBM_CUDA_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.cuh"
  )

  add_library(lib-lbm-sim-cuda INTERFACE ${LBM_CUDA_HEADERS})
  target_include_directories(lib-lbm-sim-cuda INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
  )
  target_link_libraries(lib-lbm-sim-cuda INTERFACE lib-lbm-sim)
  # target_compile_options(lib-lbm-sim-cuda INTERFACE
  #   $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
  # )
  # target_compile_definitions(lib-lbm-sim-cuda INTERFACE LBM_ENABLE_CUDA)
endif()

if (${PROJECT_IS_TOP_LEVEL})
  add_subdirectory(simulations)
endif()

add_subdirectory(docs)
