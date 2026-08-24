include(FetchContent)
FetchContent_Declare(
  quill
  GIT_REPOSITORY https://github.com/odygrd/quill
  GIT_TAG v12.1.0
)
FetchContent_MakeAvailable(quill)

set(LBM_LOG_LEVEL_VALUES
TRACE_L3 TRACE_L2 TRACE_L1 DEBUG INFO WARNING ERROR CRITICAL)

set_property(CACHE LBM_LOG_LEVEL PROPERTY STRINGS ${LBM_LOG_LEVEL_VALUES})

if(NOT DEFINED LBM_LOG_LEVEL)
  set(LBM_LOG_LEVEL_DEFAULT ${LBM_LOG_LEVEL})
endif()

if(NOT LBM_LOG_LEVEL IN_LIST LBM_LOG_LEVEL_VALUES)
  message(FATAL_ERROR
    "LBM_LOG_LEVEL must be one of:
    ${LBM_LOG_LEVEL_VALUES} (got: ${LBM_LOG_LEVEL})"
  )
endif()

set(LBM_LOG_LEVEL LBM_LOG_LEVEL_DEFAULT
CACHE STRING "Compile-time active Quill log level")

message(STATUS "LBM_LOG_LEVEL (Quill compile-time log level): ${LBM_LOG_LEVEL}")

function(lbm_apply_quill_log_level target)
  target_compile_definitions(${target} PUBLIC
    QUILL_COMPILE_ACTIVE_LOG_LEVEL=QUILL_COMPILE_ACTIVE_LOG_LEVEL_${LBM_LOG_LEVEL}
  )
  target_compile_definitions(${target} PUBLIC
    LBM_LOG_LEVEL_STR="${LBM_LOG_LEVEL}"
  )
endfunction()
