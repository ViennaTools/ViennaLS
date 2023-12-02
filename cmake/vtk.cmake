macro(setup_vtk_env TARGET OUTPUT)
  if(NOT WIN32)
    message(
      STATUS "[ViennaLS] Skipping VTK-Environment setup for ${TARGET} (Only required on Windows)")
  else()
    message(STATUS "[ViennaLS] Setting up VTK for ${TARGET}")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_BINARY_DIR}/${OUTPUT}>)

    add_custom_command(
      TARGET ${TARGET}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory ${VTK_LIBS} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif()
endmacro()
