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

function(viennals_patch_vtk_msvc_stdext VTK_SOURCE_DIR)
  if(NOT MSVC)
    return()
  endif()

  set(_vtk_fmt_header
      "${VTK_SOURCE_DIR}/ThirdParty/diy2/vtkdiy2/include/vtkdiy2/fmt/format.h")

  if(NOT EXISTS "${_vtk_fmt_header}")
    message(
      WARNING
        "[ViennaLS] Could not find VTK diy2/fmt header for MSVC stdext patch: ${_vtk_fmt_header}"
    )
    return()
  endif()

  file(READ "${_vtk_fmt_header}" _vtk_fmt_contents)

  set(_patched_guard
      "#if defined(_SECURE_SCL) && (!defined(_MSC_VER) || _MSC_VER < 1951)")

  string(FIND "${_vtk_fmt_contents}" "${_patched_guard}" _already_patched)

  if(NOT _already_patched EQUAL -1)
    message(STATUS "[ViennaLS] VTK MSVC stdext patch already applied")
    return()
  endif()

  set(_old_guard "#ifdef _SECURE_SCL")

  string(FIND "${_vtk_fmt_contents}" "${_old_guard}" _old_guard_pos)

  if(_old_guard_pos EQUAL -1)
    message(
      WARNING
        "[ViennaLS] VTK MSVC stdext patch was not applied; expected guard not found in ${_vtk_fmt_header}"
    )
    return()
  endif()

  string(REPLACE
         "${_old_guard}"
         "${_patched_guard}"
         _vtk_fmt_contents
         "${_vtk_fmt_contents}")

  file(WRITE "${_vtk_fmt_header}" "${_vtk_fmt_contents}")

  message(STATUS "[ViennaLS] Applied VTK MSVC stdext patch")
endfunction()

function(viennals_patch_vtk_openmp_nested VTK_SOURCE_DIR)
  set(_vtk_smp_openmp
      "${VTK_SOURCE_DIR}/Common/Core/SMP/OpenMP/vtkSMPTools.cxx")

  if(NOT EXISTS "${_vtk_smp_openmp}")
    message(
      WARNING
        "[ViennaLS] Could not find VTK OpenMP SMP source for omp_set_nested patch: ${_vtk_smp_openmp}"
    )
    return()
  endif()

  file(READ "${_vtk_smp_openmp}" _vtk_smp_contents)

  string(FIND "${_vtk_smp_contents}" "VIENNALS_PATCH_OMP_SET_NESTED" _already_patched)

  if(NOT _already_patched EQUAL -1)
    message(STATUS "[ViennaLS] VTK OpenMP nested-parallelism patch already applied")
    return()
  endif()

  string(FIND "${_vtk_smp_contents}" "omp_set_nested(" _has_set_nested)

  if(_has_set_nested EQUAL -1)
    message(STATUS "[ViennaLS] VTK OpenMP nested-parallelism patch not needed")
    return()
  endif()

  # Replace:
  #   omp_set_nested(isNested);
  #
  # with the modern max-active-levels control.
  #
  # Apple Clang/Homebrew libomp can expose omp_set_max_active_levels() while
  # still advertising an older _OPENMP value, so gating on the 201811 macro
  # leaves the deprecated API call in place on macOS.
  #
  # OpenMP semantics:
  #   omp_set_nested(true)  -> enable nested parallelism
  #   omp_set_nested(false) -> max-active-levels = 1
  #
  # Setting a large max-active-levels value is clamped by the runtime to the
  # supported limit, so it is a practical substitute for the "true" case.
  string(REGEX REPLACE
         "omp_set_nested\\(([^\\)]*)\\);"
         "#if defined(_OPENMP) && _OPENMP >= 200805\n"
         "  /* VIENNALS_PATCH_OMP_SET_NESTED */\n"
         "  omp_set_max_active_levels((\\1) ? 1024 : 1);\n"
         "#else\n"
         "  omp_set_nested(\\1);\n"
         "#endif"
         _vtk_smp_contents
         "${_vtk_smp_contents}")

  # Also avoid the matching deprecated query routine if VTK uses it.
  string(REGEX REPLACE
         "omp_get_nested\\(\\)"
         "/* VIENNALS_PATCH_OMP_GET_NESTED */ (omp_get_max_active_levels() > 1)"
         _vtk_smp_contents
         "${_vtk_smp_contents}")

  file(WRITE "${_vtk_smp_openmp}" "${_vtk_smp_contents}")

  message(STATUS "[ViennaLS] Applied VTK OpenMP nested-parallelism patch")
endfunction()
