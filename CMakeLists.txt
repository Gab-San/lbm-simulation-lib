cmake_minimum_required(VERSION 3.10)

project(lbm-2-lbm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE PROJECT_SOURCES "include/lbm-2-lbm/*.hpp")
file(GLOB_RECURSE PROJECT_SIMULATIONS "simulations/*.cpp")
add_executable(${PROJECT_NAME} ${PROJECT_SIMULATIONS})
#target_compile_options(
# ${PROJECT_NAME} PRIVATE "-Wnon-virtual-dtor -Wsuggest-override"
# )
target_include_directories(${PROJECT_NAME} PRIVATE "include/")
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

find_package(Doxygen)
if(Doxygen_FOUND)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs)
	message(STATUS
	"Doxygen found. Documentation will be outputted to
	${DOXYGEN_OUTPUT_DIRECTORY}")
	doxygen_add_docs(lbm_doc
	${PROJECT_SOURCES}
	ALL
	COMMENT "Generating documentation with Doxygen")
endif()
