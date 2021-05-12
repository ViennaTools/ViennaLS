## external project base directory for all projects in this hierarchy
set_property(DIRECTORY PROPERTY EP_BASE ${DEPENDENCIES_DIR})

add_subdirectory(vtk)
add_subdirectory(viennahrle)
add_subdirectory(pybind11)