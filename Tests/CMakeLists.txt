file(GLOB entries LIST_DIRECTORIES true "*")

foreach(entry ${entries})
  IF(NOT IS_DIRECTORY ${entry})
    CONTINUE()
  ENDIF()

  add_subdirectory(${entry})
endforeach()
