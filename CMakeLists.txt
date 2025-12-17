cmake_minimum_required(VERSION 3.10)

project(lbm-2-lbm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Trova OpenMP
find_package(OpenMP REQUIRED)

# Aggiungi l'eseguibile
add_executable(lbm-2-lbm
    src/main.cpp
    src/LBM.cpp
)

# Linka OpenMP al target
target_link_libraries(lbm-2-lbm PRIVATE OpenMP::OpenMP_CXX)

# Include directory per gli header
target_include_directories(lbm-2-lbm PRIVATE src)