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
  file(
    GLOB_RECURSE _vtk_smp_openmp_sources
    "${VTK_SOURCE_DIR}/Common/Core/SMP/OpenMP/*.cxx"
    "${VTK_SOURCE_DIR}/Common/Core/SMP/OpenMP/*.txx"
    "${VTK_SOURCE_DIR}/Common/Core/SMP/OpenMP/*.h")

  if(NOT _vtk_smp_openmp_sources)
    message(WARNING "[ViennaLS] Could not find VTK OpenMP SMP sources for omp_set_nested patch")
    return()
  endif()

  set(_vtk_smp_openmp_patch_count 0)

  foreach(_vtk_smp_openmp IN LISTS _vtk_smp_openmp_sources)
    file(READ "${_vtk_smp_openmp}" _vtk_smp_contents)

    string(FIND "${_vtk_smp_contents}" "omp_set_nested(" _has_set_nested)
    string(FIND "${_vtk_smp_contents}" "omp_get_nested()" _has_get_nested)

    if(_has_set_nested EQUAL -1 AND _has_get_nested EQUAL -1)
      continue()
    endif()

    if(NOT _has_set_nested EQUAL -1)
      # omp_set_nested is deprecated. Use max-active-levels directly so the
      # deprecated call is removed even when a compiler reports an old _OPENMP
      # macro while linking against a modern runtime.
      string(REGEX REPLACE
             "omp_set_nested\\(([^\\)]*)\\);"
             "/* VIENNALS_PATCH_OMP_SET_NESTED */\n"
             "  omp_set_max_active_levels((\\1) ? 1024 : 1);"
             _vtk_smp_contents
             "${_vtk_smp_contents}")
    endif()

    if(NOT _has_get_nested EQUAL -1)
      string(REGEX REPLACE
             "omp_get_nested\\(\\)"
             "/* VIENNALS_PATCH_OMP_GET_NESTED */ (omp_get_max_active_levels() > 1)"
             _vtk_smp_contents
             "${_vtk_smp_contents}")
    endif()

    file(WRITE "${_vtk_smp_openmp}" "${_vtk_smp_contents}")
    math(EXPR _vtk_smp_openmp_patch_count "${_vtk_smp_openmp_patch_count} + 1")
    message(STATUS "[ViennaLS] Patched VTK OpenMP nested-parallelism source: ${_vtk_smp_openmp}")
  endforeach()

  if(_vtk_smp_openmp_patch_count EQUAL 0)
    message(STATUS "[ViennaLS] VTK OpenMP nested-parallelism patch not needed")
  else()
    message(STATUS "[ViennaLS] Applied VTK OpenMP nested-parallelism patch")
  endif()
endfunction()
