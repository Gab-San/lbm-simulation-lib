# I main usano path relativi alla working directory: "out/..." per i frame e
# i profili, "benchmarks/ghia/..." per il confronto con Ghia et al. Vanno
# quindi lanciati da questa build dir (dove stanno anche gli eseguibili), cosi'
# l'output finisce in build/simulations/out/. Perche' anche i benchmark
# risolvano da qui, li copiamo accanto agli eseguibili in fase di configure.
if(EXISTS "${CMAKE_SOURCE_DIR}/benchmarks")
  file(COPY "${CMAKE_SOURCE_DIR}/benchmarks"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif()

file(GLOB_RECURSE CPP_SOURCES "*.cpp")

foreach(SOURCE ${CPP_SOURCES})
  get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)
  add_executable(${TARGET_NAME} ${SOURCE})
  target_link_libraries(${TARGET_NAME} PRIVATE ${LBM_SIM_LIB})
endforeach()

if(LBM_ENABLE_CUDA)
  file(GLOB_RECURSE CU_SOURCES "*.cu")
  foreach(SOURCE ${CU_SOURCES})
    get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)
    set(TARGET_NAME "cuda_${TARGET_NAME}")
    add_executable(${TARGET_NAME} ${SOURCE})
    target_link_libraries(${TARGET_NAME} PRIVATE ${LBM_SIM_LIB_CUDA})
  endforeach()
endif()
