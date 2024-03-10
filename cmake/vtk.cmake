macro(setup_vtk_env TARGET OUTPUT)
  message(STATUS "[ViennaLS] Setting up VTK-Environment for ${TARGET}")

  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory $<TARGET_FILE_DIR:vtksys> ${OUTPUT})
endmacro()
