# if this is a static build, build vtk statically
if(VIENNALS_STATIC_BUILD)
  set(VTK_EXTERNAL_BUILD_SHARED OFF)
  set(VTK_EXTERNAL_INSTALL "static-vtk-external")
  message(STATUS "Static build: Using static VTK libs.")
else()
  set(VTK_EXTERNAL_BUILD_SHARED ON)
  set(VTK_EXTERNAL_INSTALL "vtk-external")
endif()

externalproject_add(vtk-external
  GIT_REPOSITORY
    https://gitlab.kitware.com/vtk/vtk.git
  # Hash of tags/v9.0.1
  GIT_TAG #v9.0.1
    96e6fa9b3ff245e4d51d49f23d40e9ad8774e85e
  # GIT_SHALLOW 1
  TMP_DIR
    "${DEPENDENCIES_DIR}/tmp/${VTK_EXTERNAL_INSTALL}"
  STAMP_DIR
    "${DEPENDENCIES_DIR}/Stamp/${VTK_EXTERNAL_INSTALL}"
  BINARY_DIR
    "${DEPENDENCIES_DIR}/Build/${VTK_EXTERNAL_INSTALL}"
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    #-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS=${VTK_EXTERNAL_BUILD_SHARED}
    -DBUILD_TESTING=OFF
    -DBUILD_EXAMPLES=OFF
    #
    -DVTK_WRAP_PYTHON=OFF
    #-DVTK_PYTHON_VERSION=3.5
    -DVTK_LEGACY_REMOVE=ON
    -DVTK_Group_StandAlone=OFF
    -DVTK_Group_Rendering=OFF
    -DModule_vtkIOXML=ON
    -DVTK_RENDERING_BACKEND=None
  # BUILD_COMMAND
    # "cmake --build . -- -j 4"
  # No install needed
  INSTALL_COMMAND
    ""
  UPDATE_COMMAND
    ""
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
)

# Get install and build directory
externalproject_get_property(vtk-external INSTALL_DIR BINARY_DIR)

# VTK should be built
set(VIENNALS_BUILD_VTK ON CACHE BOOL "Whether VTK should be installed.")

set(
  VTK_DIR "${BINARY_DIR}"
  CACHE PATH "Path to VTK installation"
  )

# add to buildDependencies target
add_dependencies(buildDependencies vtk-external)
