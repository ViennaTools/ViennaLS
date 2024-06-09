macro(setup_vtk_env TARGET OUTPUT)
  message(STATUS "[ViennaLS] Setting up VTK-Environment for ${TARGET}")

  if(NOT TARGET vtksys)
    message(WARNING "[ViennaLS] Could not find VTK-Target")
    return()
  endif()

  # We expect all of the VTK binaries to be present in the same directory to which "vtksys" is
  # built. This is currently the case, and has been the case for prior vtk versions - However we
  # should keep an eye on this.

  add_custom_command(
    TARGET ${TARGET}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory $<TARGET_FILE_DIR:vtksys> ${OUTPUT})
endmacro()

# --------------------------------------------------------------------------------------------------------
# This is one of the most, if not the most, *cursed* workaround I've ever deployed.
# In an effort to reduce build times and to deduplicate VTK library files,
# "import_vtk_python" does the following:
#
# 1. Setup a virtual python environment and install the python vtk (=9.3.0) package to it,
#    this will also pull in pre-built vtk shared libraries as long as VTK is not installed on the system.
#
# 2. Download the VTK (=9.3.0) tar-ball from the archlinux archive project because the
#    KitWare GitLab is notriously cumbersome to download from (sometimes download speeds max out at 100KB/s)
#
# 3. Add an IMPORTED target which includes the headers from the afforementioned tarball, and links
#    against shared libraries supplied by the python package.
#
# Note: We may be able to omit this in the future when https://gitlab.kitware.com/vtk/vtk/-/issues/19279#fn-2-8090
#       is fully implemented.
# --------------------------------------------------------------------------------------------------------

macro(import_vtk_python)
  find_package(VTK QUIET)

  if(VTK_FOUND)
    message(
      FATAL_ERROR "[ViennaLS] Cannot use VTK-Python libs when system installation is present!")
  endif()

  set(VIRTUAL_ENV "${CMAKE_CURRENT_BINARY_DIR}/vtk-python-env")

  if(NOT EXISTS "${VIRTUAL_ENV}")
    message(STATUS "[ViennaLS] Setting up VTK-Python Environment")

    find_program(python-bin python REQUIRED)
    message(STATUS "[ViennaLS] Found Python: ${python-bin}")

    execute_process(COMMAND "${python-bin}" -m venv "${VIRTUAL_ENV}")
    execute_process(COMMAND "${VIRTUAL_ENV}/bin/python" -m pip install vtk==9.3.0)
  endif()

  CPMAddPackage(
    NAME VTK_TARBALL
    VERSION 9.3.0
    URL "https://archive.archlinux.org/packages/v/vtk/vtk-9.3.0-1-x86_64.pkg.tar.zst")

  if(UNIX)
    file(GLOB vtk-libraries "${VIRTUAL_ENV}/lib/python*/site-packages/vtkmodules/*.[!cpython]*.so")
  else()
    message(
      FATAL_ERROR
        "[ViennaLS] Your platform or compiler is currently not compatible with `VIENNALS_VTK_PYTHON_LIBS`"
    )
  endif()

  add_library(VTK_PythonLibs INTERFACE IMPORTED)
  add_library(VTK::PythonLibs ALIAS VTK_PythonLibs)

  target_link_libraries(VTK_PythonLibs INTERFACE ${vtk-libraries})
  target_include_directories(VTK_PythonLibs INTERFACE "${VTK_TARBALL_SOURCE_DIR}/usr/include/vtk")

  message(STATUS "[ViennaLS] Successfully created VTK::PythonLibs target")
endmacro()
