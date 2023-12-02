add_custom_target(ViennaLS_Tests ALL)

if(NOT VTK_DIR)
  include("../cmake/vtk.cmake")
  setup_vtk_env(ViennaLS_Tests "tests")
endif()

file(
  GLOB entries
  LIST_DIRECTORIES true
  "*")

foreach(entry ${entries})
  if(NOT IS_DIRECTORY ${entry})
    continue()
  endif()

  add_subdirectory(${entry})
endforeach()
