## external project base directory for all projects in this hierarchy
set_property(DIRECTORY PROPERTY EP_BASE ${DEPENDENCIES_DIR})

if(VIENNALS_USE_VTK)
    add_subdirectory(vtk)
endif()
add_subdirectory(viennahrle)
add_subdirectory(pybind11)