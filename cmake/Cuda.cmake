# lbm::sim-cuda -- the CUDA layer on top of lbm::sim.
#
# Included only from the root CMakeLists, inside if(LBM_ENABLE_CUDA): the guard
# lives with the decision, not in here.
#
# Everything this module sets globally (host compiler, architectures, separable
# compilation) must be set BEFORE enable_language(CUDA) or before the first
# CUDA target exists, which is why the order below matters. The root file
# includes this before add_subdirectory(simulations), so the cuda_* executables
# inherit it.

include_guard(GLOBAL)

set(LBM_SIM_LIB_CUDA lbm-sim-cuda)

# ---------------------------- HOST COMPILER ----------------------------------
#
# nvcc drives a host compiler for the CPU half of every .cu. Letting it pick its
# own default (usually g++) while the rest of the project is built with a
# different CXX makes a cuda_* executable mix two OpenMP runtimes: host code
# from .cu files compiled by g++ (libgomp), linked against whichever library
# OpenMP::OpenMP_CXX found for the other compiler (libomp, if CXX is clang).
# That only holds together because libomp also exports the GOMP_* symbols -- an
# accident, not a design. So align the host compiler with CMAKE_CXX_COMPILER.
#
# Must be set before enable_language(CUDA): afterwards nvcc has already been
# probed and the variable is no longer read. It stays overridable from the
# command line, for host compilers nvcc does not support.
#
# The Visual Studio generators are the exception: they always use the C++
# compiler that ships with Visual Studio, and setting this variable there only
# earns a CMake warning (see CMakeDetermineCUDACompiler.cmake). That generator
# gives MSVC to both languages anyway, so the mismatch cannot arise.
if(NOT CMAKE_GENERATOR MATCHES "Visual Studio")
  if(NOT CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
  endif()
endif()

# Also before enable_language(CUDA): from 3.18 CMake fills this in with the
# compiler's own default if it is still undefined by then, and `if(NOT ...)`
# would no longer fire.
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "Target GPU architectures")
endif()

# Relocatable device code, the safe guard for __constant__.
#
# This used to be set_target_properties(CUDA_SEPARABLE_COMPILATION ON) on an
# INTERFACE library, where it does nothing: the property is consumed when the
# *executable* is compiled and does not travel across a link edge. So the guard
# was never actually in effect. As a CMAKE_ variable it initialises the property
# on every CUDA target created after this point, including simulations/.
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

# ------------------------------- LANGUAGE ------------------------------------
#
# check_language first: enable_language(CUDA) on a machine without nvcc aborts
# with a message about a broken compiler test, which is not the useful thing to
# read when you have simply asked for a backend you cannot build.
include(CheckLanguage)
check_language(CUDA)
if(NOT CMAKE_CUDA_COMPILER)
  message(FATAL_ERROR
    "LBM_ENABLE_CUDA=ON but no CUDA compiler was found. Install the CUDA "
    "toolkit, point CMAKE_CUDA_COMPILER at nvcc, or configure with "
    "-DLBM_ENABLE_CUDA=OFF.")
endif()

enable_language(CUDA)

# -------------------------------- TARGET -------------------------------------

file(GLOB_RECURSE LBM_CUDA_HEADERS CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/include/lbm-sim/*.cuh")

add_library(${LBM_SIM_LIB_CUDA} INTERFACE)
add_library(lbm::sim-cuda ALIAS ${LBM_SIM_LIB_CUDA})

# PRIVATE on an INTERFACE library: listed for IDEs, not injected as sources into
# every consumer (which is what `add_library(... INTERFACE ${headers})` did).
target_sources(${LBM_SIM_LIB_CUDA} PRIVATE ${LBM_CUDA_HEADERS})

# Include directories, cxx_std_17, OpenMP and lbm::logging all arrive through
# this one edge. lbm-logging was previously linked three more times by hand.
target_link_libraries(${LBM_SIM_LIB_CUDA} INTERFACE lbm::sim)

target_compile_features(${LBM_SIM_LIB_CUDA} INTERFACE cuda_std_17)

# If device code needs constexpr host functions (std::array::operator[] in a
# kernel, say), this is the flag -- kept commented because it hides real errors
# when it is not needed:
# target_compile_options(${LBM_SIM_LIB_CUDA} INTERFACE
#   $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

# ---------------------------- OPENMP UNDER NVCC ------------------------------
#
# OpenMP::OpenMP_CXX exposes its flags behind $<COMPILE_LANGUAGE:CXX> (see
# FindOpenMP.cmake), so a .cu used to link the OpenMP runtime while never being
# compiled with -fopenmp: every #pragma omp in host code reachable from a .cu
# was silently ignored. Not hypothetical -- compute_boundary_mask() in
# boundaries.hpp is a host pass that runs in CUDA builds too, and it ran serial.
#
# The flag belongs to the host compiler nvcc drives, not to nvcc itself. Literal
# -fopenmp rather than ${OpenMP_CXX_FLAGS}: the latter is detected for
# CMAKE_CXX_COMPILER (aligned with the host compiler above, but overridable) and
# can name a flag the host compiler rejects -- clang gives -fopenmp=libomp,
# which g++ refuses. Plain -fopenmp is understood by both.
if(MSVC)
  set(LBM_CUDA_HOST_OPENMP_FLAG "/openmp")
else()
  set(LBM_CUDA_HOST_OPENMP_FLAG "-fopenmp")
endif()

target_compile_options(${LBM_SIM_LIB_CUDA} INTERFACE
  $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${LBM_CUDA_HOST_OPENMP_FLAG}>)

# ------------------------------- HYGIENE -------------------------------------
#
# The half-finished move: include/lbm-sim/cuda/ and include/lbm-sim/backend/cuda/
# both hold structs.cuh and utils.cuh. The glob above picks up both copies, and
# if they share include guards whichever one is reached first wins silently.
if(EXISTS "${PROJECT_SOURCE_DIR}/include/lbm-sim/cuda" AND
   EXISTS "${PROJECT_SOURCE_DIR}/include/lbm-sim/backend/cuda")
  message(WARNING
    "Both include/lbm-sim/cuda/ and include/lbm-sim/backend/cuda/ exist; the "
    "*.cuh glob is picking up both copies. Finish the move and delete one.")
endif()
