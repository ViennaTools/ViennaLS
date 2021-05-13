## external project base directory for all projects in this hierarchy
set_property(DIRECTORY PROPERTY EP_BASE ${DEPENDENCIES_DIR})

if(VIENNALS_USE_VTK)
    # if no directory is specified, or this is a re-run
    if(NOT DEFINED VTK_DIR OR VIENNALS_BUILD_VTK)
        add_subdirectory(vtk)
        message(STATUS "VTK_DIR: ${VTK_DIR}")
    endif()
endif()

# build viennahrle if no path was specified
if(NOT DEFINED ViennaHRLE_DIR OR VIENNALS_BUILD_VIENNAHRLE)
    add_subdirectory(viennahrle)
    message(STATUS "ViennaHRLE_DIR: ${ViennaHRLE_DIR}")
endif()

# build pybind11 if no path was specified
# if(NOT DEFINED pybind11_DIR OR VIENNALS_BUILD_PYBIND11)
#     add_subdirectory(pybind11)
# endif()