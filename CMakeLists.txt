cmake_minimum_required(VERSION 3.10)

project(lbm-2-lbm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp" "src/*.h" "src/*.hpp")
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})

# Trova OpenMP
find_package(OpenMP REQUIRED)

# Linka OpenMP al target
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
