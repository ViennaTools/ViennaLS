macro(setup_vtk_env TARGET OUTPUT)
  message(STATUS "[ViennaLS] Setting up VTK-Environment for ${TARGET}")

  # We expect all of the VTK binaries to be present in the same directory to which "vtksys" is
  # built. This is currently the case, and has been the case for prior vtk versions - However we
  # should keep an eye on this.

  if(TARGET vtksys)
    add_custom_command(
      TARGET ${TARGET}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory $<TARGET_FILE_DIR:vtksys> ${OUTPUT})
  else()
    add_custom_command(
      TARGET ${TARGET}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${OUTPUT})

    foreach(target ${VTK_LIBRARIES})
      get_target_property(target_path ${target} LOCATION)

      if(target_path STREQUAL "target_path-NOTFOUND")
        continue()
      endif()

      add_custom_command(
        TARGET ${TARGET}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t ${OUTPUT} ${target_path})
    endforeach()
  endif()
endmacro()
