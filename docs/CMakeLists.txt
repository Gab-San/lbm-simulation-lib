find_package(Doxygen)
if(Doxygen_FOUND)
  set(DOXYGEN_ROOT_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/doxygen")
  set(DOXYGEN_PROJECT_NAME lbm-sim)
  set(
    DOXYGEN_OUTPUT_DIRECTORY
    ${DOXYGEN_ROOT_OUTPUT_DIRECTORY}/lbm-sim
  )

  set(DOXYGEN_USE_MATHJAX YES)

  message(STATUS
    "Doxygen found. Documentation will be outputted to
    ${DOXYGEN_OUTPUT_DIRECTORY}"
  )

  set(LBM_DOXYGEN_INPUTS ${LBM_CORE_HEADERS})
  if(LBM_ENABLE_CUDA)
    list(APPEND LBM_DOXYGEN_INPUTS ${LBM_CUDA_HEADERS})
  endif()

  doxygen_add_docs(lbm-sim-docs
    ${LBM_DOXYGEN_INPUTS}
    COMMENT "Generating documentation with Doxygen"
  )


  set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
  set(DOXYGEN_OUTPUT_DIRECTORY ${DOXYGEN_ROOT_OUTPUT_DIRECTORY})
  set(
    DOXYGEN_USE_MDFILE_AS_MAINPAGE
    ${CMAKE_CURRENT_SOURCE_DIR}/mainpage.md
  )

  doxygen_add_docs(lbm-docs-full
    ${CMAKE_CURRENT_SOURCE_DIR}/mainpage.md
    COMMENT "Generating main documentation portal"
  )

  # --- THE AUTOMATION MAGIC ---
  # This makes the main target depend on the two library targets
  add_dependencies(lbm-docs-full lbm-sim-docs)

  message(STATUS
    "Doxygen automation configured. Use 'lbm-docs-full' to build all."
  )
endif()
