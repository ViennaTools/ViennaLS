# ##################################################################################################
# Check VTK Dependency
# ##################################################################################################

# Get the number of cores
include(ProcessorCount)
ProcessorCount(NPROC)

# if this is a static build, also build vtk statically
if(VIENNALS_STATIC_BUILD)
  set(VTK_EXTERNAL_BUILD_SHARED OFF)
  set(VTK_EXTERNAL_INSTALL "static_vtk_external")
  message(STATUS "Static build: Using static VTK libs.")
else()
  set(VTK_EXTERNAL_BUILD_SHARED ON)
  set(VTK_EXTERNAL_INSTALL "vtk_external")
endif()

set(VTK_COMPONENTS
    vtkCommonCore;vtkIOCore;vtkIOXML;vtkIOXMLParser;vtkFiltersCore;vtkFiltersGeneral;vtkFiltersGeometry
)

# List of VTK 9+ library targets
set(VTK_LIBRARY_TARGETS
    VTK::CommonComputationalGeometry;VTK::CommonCore;VTK::CommonDataModel;VTK::CommonExecutionModel;VTK::CommonMath;VTK::CommonMisc;VTK::CommonSystem;VTK::CommonTransforms;VTK::FiltersCore;VTK::FiltersGeneral;VTK::FiltersGeometry;VTK::IOCore;VTK::IOXML;VTK::IOXMLParser;VTK::WrappingTools
)

# If the path to an installation was provided by the user
if(DEFINED VTK_DIR AND NOT DEFINED CACHE{VTK_DIR})
  # Required for windows to convert backslashes to forward slashes
  file(TO_CMAKE_PATH "${VTK_DIR}" VTK_DIR)
  set(VTK_DIR
      "${VTK_DIR}"
      CACHE PATH "Path to VTK installation" FORCE)
  find_package(VTK CONFIG REQUIRED PATHS ${VTK_DIR} COMPONENTS "${VTK_COMPONENTS}" NO_DEFAULT_PATH)
  if(VIENNALS_STATIC_BUILD)
    check_all_type(ALL_STATIC "STATIC_LIBRARY" "${VTK_LIBRARY_TARGETS}")
    if(NOT ${ALL_STATIC})
      message(
        FATAL_ERROR
          "Not all of the targets which were exported by the provided VTK installation were built to be statically linked."
      )
    endif()
  endif()
  message(STATUS "Found VTK at ${VTK_DIR}.")
  add_library(vtk_external INTERFACE) # dummy target
else() # Automatically detect installation
  # First try to locate any local build
  find_package(VTK QUIET PATHS "${DEPENDENCIES_DIR}/Install/${VTK_EXTERNAL_INSTALL}"
               COMPONENTS "${VTK_COMPONENTS}" NO_DEFAULT_PATH)

  # If that fails look for a system installation
  if(NOT VTK_FOUND)
    find_package(VTK QUIET COMPONENTS "${VTK_COMPONENTS}")
  endif()

  if(VIENNALS_STATIC_BUILD)
    # Check if the library can be statically linked
    check_all_type(ALL_STATIC "STATIC_LIBRARY" "${VTK_LIBRARY_TARGETS}")
    if(NOT ${ALL_STATIC})
      set(VTK_FOUND FALSE)
      message(
        STATUS
          "The detected installation of VTK was not built using BUILD_SHARED_LIBS=OFF. Static VTK libraries will be built from source."
      )
    endif()
  endif()

  if(VTK_FOUND)
    message(STATUS "Found VTK at ${VTK_DIR}.")
    add_library(vtk_external INTERFACE) # dummy target
  else()
    message(
      STATUS
        "Couldn't find VTK. Dependency will be built with the buildDependencies target (e.g. `make buildDependencies`)."
    )

    ExternalProject_Add(
      vtk_external
      GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
      # Hash of tags/v9.2.6
      GIT_TAG f2c452c9c42005672a3f3ed9218dd9a7fecca79a
      GIT_SHALLOW true
      GIT_SUBMODULES ""
      TMP_DIR "${DEPENDENCIES_DIR}/tmp/${VTK_EXTERNAL_INSTALL}"
      STAMP_DIR "${DEPENDENCIES_DIR}/Stamp/${VTK_EXTERNAL_INSTALL}"
      BINARY_DIR "${DEPENDENCIES_DIR}/Build/${VTK_EXTERNAL_INSTALL}"
      INSTALL_DIR "${DEPENDENCIES_DIR}/Install/${VTK_EXTERNAL_INSTALL}"
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
                 -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                 -DCMAKE_SKIP_RPATH=OFF
                 -DCMAKE_SKIP_INSTALL_RPATH=ON
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
                 -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
                 -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
                 -DBUILD_SHARED_LIBS:BOOL=${VTK_EXTERNAL_BUILD_SHARED}
                 -DVTK_BUILD_TESTING=OFF
                 -DVTK_BUILD_EXAMPLES=OFF
                 -DVTK_BUILD_DOCUMENTATION=OFF
                 -DVTK_VERSIONED_INSTALL=OFF
                 -DVTK_WRAP_PYTHON:BOOL=NO
                 -DVTK_SMP_IMPLEMENTATION_TYPE:STRING="OpenMP"
                 -DVTK_GROUP_ENABLE_Rendering:BOOL=NO
                 -DVTK_GROUP_ENABLE_Imaging:BOOL=NO
                 -DVTK_GROUP_ENABLE_MPI:BOOL=NO
                 -DVTK_GROUP_ENABLE_Qt:BOOL=NO
                 -DVTK_GROUP_ENABLE_StandAlone:BOOL=NO
                 -DVTK_GROUP_ENABLE_Views:BOOL=NO
                 -DVTK_GROUP_ENABLE_Web:BOOL=NO
                 -DVTK_LEGACY_REMOVE:BOOL=ON
                 -DVTK_MODULE_ENABLE_VTK_CommonExecutionModel:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonMisc:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonSystem:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonMath:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonCore:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonTransforms:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_CommonComputationalGeometry:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_IOCore:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_IOXMLParser:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_IOXML:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_FiltersCore:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_FiltersGeneral:BOOL=YES
                 -DVTK_MODULE_ENABLE_VTK_FiltersGeometry:BOOL=YES
                 -Wno-dev
      BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NPROC} --config Release
      INSTALL_COMMAND ${CMAKE_COMMAND} --install .
      UPDATE_COMMAND ""
      USES_TERMINAL_DOWNLOAD 1
      USES_TERMINAL_UPDATE 1
      GIT_PROGRESS 1)

    # Get install directory
    ExternalProject_Get_Property(vtk_external INSTALL_DIR)

    # Set the required variables
    set(VTK_DIR "${INSTALL_DIR}")

    set(VTK_VERSION 9.2)
    set(VTK_INCLUDE_DIRS "${VTK_DIR}/include/vtk-${VTK_VERSION}")

    # add to buildDependencies target
    add_dependencies(buildDependencies vtk_external)
  endif()

  # Write VTK_DIR, VTK_VERSION and VTK_INCLUDE_DIRS to the cache so that we can use it in the
  # toplevel CMakeLists file
  set(VTK_DIR
      ${VTK_DIR}
      CACHE PATH "Path to VTK installation" FORCE)

  set(VTK_VERSION
      ${VTK_VERSION}
      CACHE STRING "The VTK version" FORCE)

  set(VTK_INCLUDE_DIRS
      ${VTK_INCLUDE_DIRS}
      CACHE STRING "Path to VTK headers" FORCE)
endif()
