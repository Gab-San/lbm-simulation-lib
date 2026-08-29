# I main usano path relativi alla working directory: "out/..." per i frame e
# i profili, "benchmarks/ghia/..." per il confronto con Ghia et al. Vanno
# quindi lanciati da questa build dir (dove stanno anche gli eseguibili), cosi'
# l'output finisce in build/simulations/out/. Perche' anche i benchmark
# risolvano da qui, li copiamo accanto agli eseguibili in fase di configure.
#
# PROJECT_SOURCE_DIR e non CMAKE_SOURCE_DIR: il secondo punta alla radice
# dell'albero CMake, che non e' questo progetto quando lo si add_subdirectory().
if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks")
  file(COPY "${PROJECT_SOURCE_DIR}/benchmarks"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# I dati di input (configurazioni .toml e benchmark di Ghia) stanno nei
# sorgenti, non nella build. Passiamo i loro path assoluti come define, cosi'
# gli eseguibili li trovano da qualunque working directory. Non li copiamo:
# puntando alla cartella sorgente, modificare un .toml ha effetto subito,
# senza ricompilare ne' rilanciare cmake.
function(lbm_add_data_paths TARGET_NAME)
  target_compile_definitions(${TARGET_NAME} PRIVATE
    LBM_CONFIGS_DIR="${PROJECT_SOURCE_DIR}/configs"
    LBM_BENCHMARKS_DIR="${PROJECT_SOURCE_DIR}/benchmarks")
endfunction()

function(lbm_add_simulation SOURCE PREFIX)
  get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)
  set(TARGET_NAME "${PREFIX}${TARGET_NAME}")

  # Il glob e' ricorsivo: due sorgenti omonime in cartelle diverse darebbero un
  # errore di target duplicato molto piu' avanti, e non su questa riga.
  if(TARGET ${TARGET_NAME})
    message(FATAL_ERROR
      "Nome di target duplicato '${TARGET_NAME}' (da ${SOURCE}): due sorgenti "
      "sotto simulations/ hanno lo stesso nome file.")
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
