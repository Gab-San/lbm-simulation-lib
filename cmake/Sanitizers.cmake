# Sanitizer support.
#
# Configure with -DLBM_SANITIZE=address (or "address;undefined", or
# "address,undefined" -- both spellings work).
#
# Why this is a global, directory-scope thing and not an INTERFACE target:
# the MSVC STL emits `#pragma detect_mismatch` directives (annotate_string,
# annotate_vector, annotate_optional) so the linker rejects a binary that mixes
# ASan-instrumented and uninstrumented translation units -- it is a real ODR
# violation, since std::vector<T>::push_back is a different function in each.
# lbm-logging and quill are static libraries fetched/built inside this project,
# so they must carry the same flags as the executables that link them. The only
# way to reach FetchContent'd subprojects is add_compile_options() before the
# corresponding add_subdirectory(), which is why the root CMakeLists includes
# this module *before* include(Logging) and before the toml++ fetch.
#
# Outputs (read by the root CMakeLists):
#   LBM_SANITIZE_ACTIVE           TRUE when any sanitizer is on
#   LBM_SANITIZE_DISABLES_OPENMP  TRUE when OpenMP must be dropped for this build

include_guard(GLOBAL)

set(LBM_SANITIZE "" CACHE STRING
    "Sanitizers to build with: address, undefined, thread, leak. \
Comma- or semicolon-separated; empty disables sanitizers.")

set(LBM_SANITIZE_ACTIVE          FALSE)
set(LBM_SANITIZE_DISABLES_OPENMP FALSE)

if(NOT LBM_SANITIZE)
  return()
endif()

# ------------------------------ NORMALISE ------------------------------------

string(REPLACE "," ";" _lbm_san_list "${LBM_SANITIZE}")
list(TRANSFORM _lbm_san_list STRIP)
list(TRANSFORM _lbm_san_list TOLOWER)
list(REMOVE_ITEM _lbm_san_list "")
list(REMOVE_DUPLICATES _lbm_san_list)

set(_lbm_san_known address undefined thread leak memory)
foreach(_san IN LISTS _lbm_san_list)
  if(NOT _san IN_LIST _lbm_san_known)
    message(FATAL_ERROR
      "LBM_SANITIZE: unknown sanitizer '${_san}'. "
      "Known values: ${_lbm_san_known}.")
  endif()
endforeach()

# ThreadSanitizer maps the address space differently and cannot coexist with
# the address/leak shadow. MemorySanitizer is Clang-only and needs every
# dependency (including libstdc++) rebuilt under it, which FetchContent will
# not give you here -- refuse it rather than produce a wall of false positives.
if("thread" IN_LIST _lbm_san_list)
  foreach(_bad address leak memory)
    if(_bad IN_LIST _lbm_san_list)
      message(FATAL_ERROR
        "LBM_SANITIZE: 'thread' cannot be combined with '${_bad}'. "
        "Configure a separate build directory for each.")
    endif()
  endforeach()
endif()

set(LBM_SANITIZE_ACTIVE TRUE)

# ================================== MSVC =====================================

if(MSVC)

  if(NOT _lbm_san_list STREQUAL "address")
    message(FATAL_ERROR
      "LBM_SANITIZE=${LBM_SANITIZE}: MSVC only implements the address "
      "sanitizer. Use -DLBM_SANITIZE=address here, and run the other "
      "sanitizers under GCC or Clang.")
  endif()

  # MSVC ASan is documented as incompatible with /openmp. Leaving the flag on
  # makes cl.exe fail with MSB8059 / an equivalent command-line error, so drop
  # OpenMP for this build and let the caller know the CPU solvers are serial.
  set(LBM_SANITIZE_DISABLES_OPENMP TRUE)

  # nvcc does not forward /fsanitize=address to the host compiler, so every .cu
  # would come out uninstrumented and collide with lbm-logging at link time --
  # that is exactly the annotate_string/annotate_vector LNK2038 wall. There is
  # no supported way to instrument nvcc host code on MSVC, so refuse the
  # combination instead of producing a build that cannot link.
  if(LBM_ENABLE_CUDA)
    message(FATAL_ERROR
      "LBM_SANITIZE=address and LBM_ENABLE_CUDA=ON are not buildable together "
      "with MSVC: nvcc will not pass /fsanitize=address to the host compiler, "
      "and mixing instrumented and uninstrumented objects fails at link with "
      "LNK2038 'annotate_string' / 'annotate_vector' / 'annotate_optional'.\n"
      "  Reconfigure with -DLBM_ENABLE_CUDA=OFF to sanitize the CPU solvers, "
      "or sanitize the CUDA build on Linux, where the host flags can be "
      "forwarded with -Xcompiler.")
  endif()

  # /Zi + /DEBUG: without debug info ASan reports raw addresses instead of a
  # call stack (that is warning C5072 / LNK4302).
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:/fsanitize=address>"
    "$<$<COMPILE_LANGUAGE:CXX>:/Zi>")

  add_link_options(/DEBUG /INCREMENTAL:NO)

  # /RTC is incompatible with ASan and CMake puts /RTC1 in the Debug flags by
  # default. Strip it from every configuration that might carry it.
  foreach(_cfg "" _DEBUG _RELWITHDEBINFO)
    string(REGEX REPLACE "/RTC(su|[1su])( |$)" ""
           CMAKE_CXX_FLAGS${_cfg} "${CMAKE_CXX_FLAGS${_cfg}}")
    string(REGEX REPLACE "/INCREMENTAL( |$)" ""
           CMAKE_EXE_LINKER_FLAGS${_cfg} "${CMAKE_EXE_LINKER_FLAGS${_cfg}}")
  endforeach()

  # Do not fold this into a target: FetchContent'd subprojects (quill) never
  # see target properties, and they must be instrumented too.
  message(STATUS
    "Sanitizers: MSVC AddressSanitizer enabled (OpenMP and CUDA disabled).")
  message(STATUS
    "  the ASan runtime DLL (clang_rt.asan_dynamic-x86_64.dll) ships with "
    "Visual Studio; run from a Developer Prompt if the exe fails to start.")

  return()
endif()

# ============================== GCC / CLANG ==================================

# One -fsanitize= flag per sanitizer rather than a comma-separated list: a
# comma inside a $<...> generator expression is parsed as an argument
# separator, and these flags have to survive being wrapped for nvcc below.
set(_lbm_san_flags "")
foreach(_san IN LISTS _lbm_san_list)
  list(APPEND _lbm_san_flags "-fsanitize=${_san}")
endforeach()

# -g so the reports carry line numbers; -fno-omit-frame-pointer so the stacks
# are usable at -O2, which is where you will actually reproduce this.
list(APPEND _lbm_san_flags -fno-omit-frame-pointer -g)

foreach(_flag IN LISTS _lbm_san_flags)
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${_flag}>")
endforeach()

if(LBM_ENABLE_CUDA)
  # The flags belong to the host compiler nvcc drives, not to nvcc. Same
  # reasoning as LBM_CUDA_HOST_OPENMP_FLAG in Cuda.cmake -- and the same
  # requirement: host code inside .cu must be instrumented identically to the
  # rest, or the CUDA executables link against a differently-built lbm-logging.
  foreach(_flag IN LISTS _lbm_san_flags)
    add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_flag}>")
  endforeach()
endif()

# Link options are per link language: a CUDA executable is linked by nvcc,
# which needs the flag wrapped, while a plain C++ executable is linked by the
# host compiler directly.
foreach(_san IN LISTS _lbm_san_list)
  add_link_options(
    "$<$<LINK_LANGUAGE:CXX>:-fsanitize=${_san}>"
    "$<$<LINK_LANGUAGE:CUDA>:-Xcompiler=-fsanitize=${_san}>")
endforeach()

message(STATUS "Sanitizers: ${_lbm_san_list}")

if(LBM_ENABLE_CUDA AND "address" IN_LIST _lbm_san_list)
  message(STATUS
    "  CUDA + ASan needs ASAN_OPTIONS=protect_shadow_gap=0 at runtime, "
    "otherwise CUDA context creation fails.")
  message(STATUS
    "  LeakSanitizer reports the CUDA runtime's own allocations; add "
    "detect_leaks=0 to ASAN_OPTIONS or use a suppression file.")
endif()
