# external project base directory for all projects in this hierarchy. Set this variable here so that
# pythonBindings are not affected
set_property(DIRECTORY PROPERTY EP_BASE ${DEPENDENCIES_DIR})

# Macro to check if all of the libraries provided by the target list LIBRARY_TARGETS are compiled as
# the provided type (TARGET_TYPE: e.g. "STATIC_LIBRARY", "SHARED_LIBRARY", ...)
macro(CHECK_ALL_TYPE RESULT TARGET_TYPE LIBRARY_TARGETS)
  set(_all_match TRUE)
  foreach(_lib ${LIBRARY_TARGETS})
    if(TARGET ${_lib})
      get_target_property(_target_type ${_lib} TYPE)
      # message(STATUS "${_lib}: ${_target_type} ${TARGET_TYPE}")
      if(NOT _target_type STREQUAL ${TARGET_TYPE})
        set(_all_match FALSE)
        break()
      endif()
    else()
      message(WARNING "CHECK_ALL_TYPE: Not all targets were found!\- caused by ${_lib}")
      set(_all_match FALSE)
      break()
    endif()
  endforeach()
  set(${RESULT} ${_all_match})
endmacro()

# OpenMP Support is required
find_package(OpenMP REQUIRED)

# ViennaHRLE is also mandatory
add_subdirectory(viennahrle)

if(VIENNALS_USE_VTK)
  add_subdirectory(vtk)
endif()

if(VIENNALS_BUILD_PYTHON)
  add_subdirectory(pybind11)
endif()
