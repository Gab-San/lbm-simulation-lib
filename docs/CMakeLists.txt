# Doxygen portal for lbm-simulation-lib.
#
# Two targets:
#   lbm-sim-docs  -- the library headers alone, no narrative pages. Useful
#                    when only the API reference is wanted.
#   lbm-docs-full -- the whole portal: mainpage, the topic pages under pages/
#                    and the same headers.
#
# doxygen_add_docs() snapshots the DOXYGEN_* variables at call time, so each
# block below sets what its own target needs *before* creating it.

find_package(Doxygen)

if(NOT Doxygen_FOUND)
  message(STATUS
    "Doxygen not found: lbm-sim-docs / lbm-docs-full are not configured.")
  return()
endif()

set(DOXYGEN_ROOT_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/doxygen")

# ------------------------------- INPUTS -------------------------------------
#
# LBM_CORE_HEADERS comes from the root CMakeLists and LBM_CUDA_HEADERS from
# cmake/Cuda.cmake, which only runs with LBM_ENABLE_CUDA=ON. Neither glob
# covers include/lbm/ (the format and logging helpers shared with the rest of
# the project), so it is globbed here.
file(GLOB_RECURSE LBM_DOXYGEN_SHARED_HEADERS CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/include/lbm/*.hpp")

set(LBM_DOXYGEN_HEADERS
  ${LBM_CORE_HEADERS}
  ${LBM_DOXYGEN_SHARED_HEADERS})

if(LBM_ENABLE_CUDA)
  list(APPEND LBM_DOXYGEN_HEADERS ${LBM_CUDA_HEADERS})
endif()

# Narrative pages, in the order they should appear in the tree. mainpage.md
# comes first: DOXYGEN_USE_MDFILE_AS_MAINPAGE only has an effect on a file that
# is also part of the input list.
set(LBM_DOXYGEN_PAGES
  "${CMAKE_CURRENT_SOURCE_DIR}/mainpage.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/architecture.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/configuration.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/output-formats.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/validation.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/performance.md"
  "${CMAKE_CURRENT_SOURCE_DIR}/pages/extending.md")

# ---------------------------- COMMON SETTINGS -------------------------------

# The library is header-only and heavily templated: the interesting entities
# are classes and function templates in headers, and there are no .cpp files
# to hide.
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_MARKDOWN_SUPPORT YES)
set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
set(DOXYGEN_QT_AUTOBRIEF YES)
set(DOXYGEN_SORT_MEMBER_DOCS NO)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)

# Every header carries an @file block and the public entities are documented,
# so the reference shows what was actually written rather than padding itself
# with bare signatures. Set WARN_IF_UNDOCUMENTED=YES below to find the gaps
# that remain.
set(DOXYGEN_EXTRACT_ALL NO)
set(DOXYGEN_EXTRACT_STATIC YES)

# Implementation details that only add noise to the reference.
set(DOXYGEN_EXCLUDE_SYMBOLS "detail::*" "assertion::*")

# LBM_HD_FUNC, RESTRICT and UNROLL_FULL are compiler/CUDA decorations: expanded
# to nothing they leave a signature Doxygen can parse and a reader can read.
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
set(DOXYGEN_PREDEFINED
  "LBM_HD_FUNC="
  "RESTRICT="
  "UNROLL_FULL="
  "__host__="
  "__device__="
  "__constant__="
  "__global__=")

# Figures referenced from the pages: the older benchmark plots next to the
# README, and the result figures the scaling and validation pages point at.
set(DOXYGEN_IMAGE_PATH
  "${PROJECT_SOURCE_DIR}/images"
  "${CMAKE_CURRENT_SOURCE_DIR}/imgResults")

# The interactive velocity-set viewers are copied verbatim into the HTML
# output so the pages can link them with a relative path.
set(DOXYGEN_HTML_EXTRA_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/assets/lbm_d3q19_directions.html"
  "${CMAKE_CURRENT_SOURCE_DIR}/assets/lbm_d3q27_directions.html")

# ------------------------------- TARGETS ------------------------------------

set(DOXYGEN_PROJECT_NAME lbm-sim)
set(DOXYGEN_PROJECT_BRIEF "Lattice Boltzmann simulation core")
set(DOXYGEN_OUTPUT_DIRECTORY "${DOXYGEN_ROOT_OUTPUT_DIRECTORY}/lbm-sim")

doxygen_add_docs(lbm-sim-docs
  ${LBM_DOXYGEN_HEADERS}
  COMMENT "Generating the lbm-sim API reference with Doxygen")

set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
set(DOXYGEN_PROJECT_BRIEF "Header-only lattice Boltzmann framework (OpenMP / CUDA)")
set(DOXYGEN_OUTPUT_DIRECTORY "${DOXYGEN_ROOT_OUTPUT_DIRECTORY}")
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE
  "${CMAKE_CURRENT_SOURCE_DIR}/mainpage.md")

doxygen_add_docs(lbm-docs-full
  ${LBM_DOXYGEN_PAGES}
  ${LBM_DOXYGEN_HEADERS}
  COMMENT "Generating the full documentation portal with Doxygen")

message(STATUS
  "Doxygen found. Documentation targets: lbm-sim-docs, lbm-docs-full "
  "(output under ${DOXYGEN_ROOT_OUTPUT_DIRECTORY})")
