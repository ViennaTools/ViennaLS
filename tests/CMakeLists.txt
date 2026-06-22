add_custom_target(ViennaLS_Tests ALL)

if(WIN32)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_BINARY_DIR}/tests>)
  viennacore_setup_vtk_env(ViennaLS_Tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

viennacore_add_subdirs(${CMAKE_CURRENT_LIST_DIR})

if(VIENNALS_ENABLE_SANITIZER)
  get_directory_property(_test_subdirs SUBDIRECTORIES)
  foreach(_subdir IN LISTS _test_subdirs)
    get_property(_subdir_tests DIRECTORY "${_subdir}" PROPERTY TESTS)
    if(_subdir_tests)
      set_tests_properties(${_subdir_tests} PROPERTIES
        ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/lsan_suppressions.txt")
    endif()
  endforeach()
endif()
