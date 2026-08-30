include_guard(GLOBAL)

include(FetchContent)
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
  )
FetchContent_MakeAvailable(tomlplusplus)

set(LBM_CONFIG_LIB lbm-config)

file(GLOB_RECURSE LBM_CONFIG_HEADERS CONFIGURE_DEPENDS
"${CMAKE_CURRENT_SOURCE_DIR}/include/lbm/config/*.hpp"
)

file(GLOB_RECURSE LBM_CONFIG_SRC CONFIGURE_DEPENDS
"${CMAKE_CURRENT_SOURCE_DIR}/src/lbm/config/*.cpp"
)

add_library(${LBM_CONFIG_LIB} STATIC ${LBM_CONFIG_HEADERS} ${LBM_CONFIG_SRC})
add_library(lbm::config ALIAS ${LBM_CONFIG_LIB})

target_include_directories(${LBM_CONFIG_LIB} PUBLIC
                           "${PROJECT_SOURCE_DIR}/include")

target_compile_features(${LBM_CONFIG_LIB} PUBLIC cxx_std_17)

target_link_libraries(${LBM_CONFIG_LIB} PRIVATE tomlplusplus::tomlplusplus)

# The rest of the project is header-only, so this is the one archive everything
# links. Keep it usable from a shared consumer.
set_target_properties(${LBM_CONFIG_LIB} PROPERTIES
                      POSITION_INDEPENDENT_CODE ON)
