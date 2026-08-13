cmake_minimum_required(VERSION 3.10)

project(lbm-simulation-lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

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

include(FetchContent)
FetchContent_Declare(
  quill
  GIT_REPOSITORY https://github.com/odygrd/quill
  GIT_TAG v12.1.0
)
FetchContent_MakeAvailable(quill)

set(LBM_LOGGING_HELPER_LIB lbm-logging)

file(
  GLOB_RECURSE
  LBM_LOGGING_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/lbm/*.cpp"
)

add_library(${LBM_LOGGING_HELPER_LIB} STATIC ${LBM_LOGGING_FILES})
target_include_directories(
  ${LBM_LOGGING_HELPER_LIB}
  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(${LBM_LOGGING_HELPER_LIB} PUBLIC quill::quill)

target_link_libraries(lib-lbm-sim INTERFACE ${LBM_LOGGING_HELPER_LIB})

if (${PROJECT_IS_TOP_LEVEL})
  add_subdirectory(simulations)
endif()

find_package(Doxygen)
if(Doxygen_FOUND)
  set(DOXYGEN_ROOT_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs/doxygen")
  set(DOXYGEN_PROJECT_NAME lbm-sim)
  set(
    DOXYGEN_OUTPUT_DIRECTORY
    ${DOXYGEN_ROOT_OUTPUT_DIRECTORY}/lbm-sim
  )

  message(STATUS
    "Doxygen found. Documentation will be outputted to
    ${DOXYGEN_OUTPUT_DIRECTORY}"
  )

  doxygen_add_docs(lbm-sim-docs
    ${LBM_CORE_HEADERS}
    COMMENT "Generating documentation with Doxygen"
  )

  set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
  set(DOXYGEN_OUTPUT_DIRECTORY ${DOXYGEN_ROOT_OUTPUT_DIRECTORY})
  set(
    DOXYGEN_USE_MDFILE_AS_MAINPAGE
    ${CMAKE_CURRENT_SOURCE_DIR}/docs/mainpage.md
  )

  doxygen_add_docs(lbm-docs-full
    ${CMAKE_CURRENT_SOURCE_DIR}/docs/mainpage.md
    COMMENT "Generating main documentation portal"
  )

  # --- THE AUTOMATION MAGIC ---
  # This makes the main target depend on the two library targets
  add_dependencies(lbm-docs-full lbm-sim-docs)

  message(STATUS
    "Doxygen automation configured. Use 'lbm-docs-full' to build all."
  )
endif()
