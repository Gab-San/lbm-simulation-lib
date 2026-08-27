if(LBM_ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED True)
  set(LBM_SIM_LIB_CUDA lbm-sim-cuda)

  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)   # adjust to your target GPU(s)
  endif()

  file(GLOB_RECURSE LBM_CUDA_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/lbm-sim/*.cuh"
  )

  add_library(${LBM_SIM_LIB_CUDA} INTERFACE ${LBM_CUDA_HEADERS})
  target_include_directories(${LBM_SIM_LIB_CUDA} INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
  )
  target_link_libraries(${LBM_SIM_LIB_CUDA} INTERFACE lbm-sim)
  # target_compile_options(lib-lbm-sim-cuda INTERFACE
  #   $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
  # )
  # target_compile_definitions(lib-lbm-sim-cuda INTERFACE LBM_ENABLE_CUDA)
  target_link_libraries(${LBM_SIM_LIB_CUDA} INTERFACE lbm-logging)

  # NOTE: Suppressing warnings on quill
  target_compile_options(${LBM_SIM_LIB_CUDA} INTERFACE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe --diag_suppress=128,--diag_suppress=2417>
  )
  target_link_libraries(${LBM_SIM_LIB_CUDA} INTERFACE ${LBM_LOGGING_HELPER_LIB})

  # Enables relocatable device code (safe guard for __constant__)
  set_target_properties(${LBM_SIM_LIB_CUDA}
  PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

