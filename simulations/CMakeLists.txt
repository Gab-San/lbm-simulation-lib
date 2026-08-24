# Auto-compile all generated simulation files in this directory
file(GLOB SOURCES "*.cpp")

foreach(SOURCE ${SOURCES})
  get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)
  add_executable(${TARGET_NAME} ${SOURCE})
  target_link_libraries(${TARGET_NAME} PRIVATE ${LBM_SIM_LIB})
endforeach()
