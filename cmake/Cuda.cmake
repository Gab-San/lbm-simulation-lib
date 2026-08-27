if(LBM_ENABLE_CUDA)
  # nvcc drives a host compiler for the CPU half of every .cu. Letting it pick
  # its own default (usually g++) while the rest of the project is built with a
  # different CXX makes a cuda_* executable mix two OpenMP runtimes: host code
  # from .cu files compiled by g++ (libgomp), linked against whichever library
  # OpenMP::OpenMP_CXX found for the other compiler (libomp, if CXX is clang).
  # That only holds together because libomp also exports the GOMP_* symbols --
  # an accident, not a design. So align the host compiler with CMAKE_CXX_COMPILER.
  #
  # This must be set BEFORE enable_language(CUDA): afterwards nvcc has already
  # been probed and the variable is no longer read. It stays overridable from
  # the command line, for host compilers nvcc does not support.
  #
  # The Visual Studio generators are the exception: they always use the C++
  # compiler that comes with Visual Studio, and setting this variable there
  # only earns a CMake warning (see CMakeDetermineCUDACompiler.cmake). Since
  # that generator gives MSVC to both languages anyway, the mismatch this
  # guards against cannot arise, so skip it.
  if(NOT CMAKE_GENERATOR MATCHES "Visual Studio")
    if(NOT CMAKE_CUDA_HOST_COMPILER)
      set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
    endif()
  endif()

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

  # OpenMP in nvcc's host pass.
  #
  # OpenMP::OpenMP_CXX exposes its flags behind $<COMPILE_LANGUAGE:CXX> (see
  # FindOpenMP.cmake), so a .cu used to link the OpenMP runtime while never
  # being compiled with -fopenmp: every #pragma omp in host code reachable from
  # a .cu was silently ignored. That is not hypothetical -- the boundary mask
  # (compute_boundary_mask in boundaries.hpp) is a host pass that runs in CUDA
  # builds too, and it ran serial.
  #
  # The flag belongs to the host compiler nvcc drives, not to nvcc itself.
  # We use a literal -fopenmp rather than ${OpenMP_CXX_FLAGS}: the latter is
  # detected for CMAKE_CXX_COMPILER, which we align with the host compiler
  # above but which stays overridable, and it can name a flag the host compiler
  # rejects (clang gives -fopenmp=libomp, which g++ refuses). Plain -fopenmp is
  # understood by both gcc and clang.
  if(MSVC)
    set(LBM_CUDA_HOST_OPENMP_FLAG "/openmp:experimental")
  else()
    set(LBM_CUDA_HOST_OPENMP_FLAG "-fopenmp")
  endif()

  target_compile_options(${LBM_SIM_LIB_CUDA} INTERFACE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${LBM_CUDA_HOST_OPENMP_FLAG}>
  )
  target_link_libraries(${LBM_SIM_LIB_CUDA} INTERFACE ${LBM_LOGGING_HELPER_LIB})

  # Enables relocatable device code (safe guard for __constant__)
  set_target_properties(${LBM_SIM_LIB_CUDA}
  PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
