# lbm::config -- TOML configuration parsing (toml++).
#
# Replaces cmake/Toml.cmake, which only did the FetchContent and left toml++ on
# the include path of everything that touched the library.
#
# The target always exists, so consumers can link it unconditionally. With
# LBM_ENABLE_CONFIG=OFF it is an empty INTERFACE library: nothing is fetched,
# include/lbm-sim/config/ is not part of the build, and anything that #includes
# those headers is expected not to be built either.

include_guard(GLOBAL)

set(LBM_CONFIG_LIB lbm-config)

add_library(${LBM_CONFIG_LIB} INTERFACE)
add_library(lbm::config ALIAS ${LBM_CONFIG_LIB})

if(LBM_ENABLE_CONFIG)
  include(FetchContent)
  FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(tomlplusplus)

  file(GLOB_RECURSE LBM_CONFIG_HEADERS CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/lbm-sim/config/*.hpp")

  target_sources(${LBM_CONFIG_LIB} PRIVATE ${LBM_CONFIG_HEADERS})

  target_include_directories(${LBM_CONFIG_LIB} INTERFACE
    "${PROJECT_SOURCE_DIR}/include")

  target_compile_features(${LBM_CONFIG_LIB} INTERFACE cxx_std_17)

  target_link_libraries(${LBM_CONFIG_LIB} INTERFACE
    tomlplusplus::tomlplusplus)

  # For code that can do without a config file: #if LBM_ENABLE_CONFIG.
  target_compile_definitions(${LBM_CONFIG_LIB} INTERFACE LBM_ENABLE_CONFIG)
endif()

message(STATUS "LBM config: ${LBM_ENABLE_CONFIG}")

# Deliberately NOT linked into lbm::sim. The core is a header-only library that
# mostly does not parse anything; toml++ was on the include path of every one of
# its consumers. The mains link lbm::config themselves, which is what makes
# -DLBM_ENABLE_CONFIG=OFF a build you can actually take.
#
# To go back to the old behaviour, one line in the root CMakeLists:
#   target_link_libraries(${LBM_SIM_LIB} INTERFACE lbm::config)
