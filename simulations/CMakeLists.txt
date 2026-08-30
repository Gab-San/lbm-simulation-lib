# The mains use paths relative to the working directory: "out/..." for frames
# and profiles, "benchmarks/ghia/..." for the comparison against Ghia et al.
# They therefore have to be launched from this build dir (which also holds the
# executables), so that the output lands in build/simulations/out/. So that
# the benchmarks resolve from here too, we copy them next to the executables
# at configure time.
#
# PROJECT_SOURCE_DIR and not CMAKE_SOURCE_DIR: the latter points to the root
# of the CMake tree, which is not this project when it is add_subdirectory()'d.
if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks")
  file(COPY "${PROJECT_SOURCE_DIR}/benchmarks"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# The input data (.toml configurations and Ghia benchmarks) lives in the
# sources, not in the build. We pass their absolute paths as defines, so the
# executables find them from any working directory. We do not copy them:
# pointing at the source directory means editing a .toml takes effect right
# away, without recompiling or re-running cmake.
function(lbm_add_data_paths TARGET_NAME)
  target_compile_definitions(${TARGET_NAME} PRIVATE
    LBM_CONFIGS_DIR="${PROJECT_SOURCE_DIR}/configs"
    LBM_BENCHMARKS_DIR="${PROJECT_SOURCE_DIR}/benchmarks")
endfunction()

function(lbm_add_simulation SOURCE PREFIX)
  get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)
  set(TARGET_NAME "${PREFIX}${TARGET_NAME}")

  # The glob is recursive: two same-named sources in different directories
  # would give a duplicate-target error much later, and not on this line.
  if(TARGET ${TARGET_NAME})
    message(FATAL_ERROR
      "Duplicate target name '${TARGET_NAME}' (from ${SOURCE}): two sources "
      "under simulations/ have the same file name.")
  endif()

  add_executable(${TARGET_NAME} ${SOURCE})
  target_link_libraries(${TARGET_NAME} PRIVATE ${ARGN} lbm::sim)
  lbm_add_data_paths(${TARGET_NAME})
endfunction()

file(GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS "*.cpp")

foreach(SOURCE ${CPP_SOURCES})
  lbm_add_simulation("${SOURCE}" "" lbm::sim)
endforeach()

if(LBM_ENABLE_CUDA)
  file(GLOB_RECURSE CU_SOURCES CONFIGURE_DEPENDS "*.cu")
  foreach(SOURCE ${CU_SOURCES})
    lbm_add_simulation("${SOURCE}" "cuda_" lbm::sim-cuda)
  endforeach()
endif()
