add_custom_target(ViennaLS_Examples ALL)

if(WIN32)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_BINARY_DIR}/examples>)
  viennacore_setup_vtk_env(ViennaLS_Examples ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

viennacore_add_subdirs(${CMAKE_CURRENT_LIST_DIR})
