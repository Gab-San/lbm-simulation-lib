cmake_minimum_required(VERSION 3.10)

project(lbm-2-lbm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Aggiungi l'eseguibile
add_executable(lbm-2-lbm
    src/main.cpp
    src/LBM.cpp
)

# Include directory per l'header dentro src/
target_include_directories(lbm-2-lbm PRIVATE src)