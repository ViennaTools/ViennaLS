add_custom_target(ViennaLS_Tests ALL)

if(WIN32 AND VIENNALS_VTK_ENV)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_BINARY_DIR}/tests>)
  setup_vtk_env(ViennaLS_Tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

viennacore_add_subdirs(${CMAKE_CURRENT_LIST_DIR})
