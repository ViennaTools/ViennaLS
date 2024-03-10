add_custom_target(ViennaLS_Tests ALL)

if(WIN32 AND NOT VIENNALS_SYSTEM_VTK)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_BINARY_DIR}/tests>)
  setup_vtk_env(ViennaLS_Tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

file(
  GLOB entries
  LIST_DIRECTORIES true
  "*")

foreach(entry ${entries})
  if(NOT IS_DIRECTORY ${entry})
    continue()
  endif()

  message(STATUS "[ViennaLS] Adding Test ${entry}")
  add_subdirectory(${entry})
endforeach()
