cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(ViennaLS VERSION 2.0.0)

set(VIENNALS_VERSION
    ${PROJECT_VERSION}
    CACHE STRING "The version of ViennaLS" FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# c++17 for inline constexpr variables
set(CMAKE_CXX_STANDARD "17")

add_definitions(-DVIENNALS_VERSION=${PROJECT_VERSION})

# set default build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

# tell VS to export all symbols to its dll files
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS
      TRUE
      CACHE BOOL "Export all symbols" FORCE)
  add_compile_definitions(_USE_MATH_DEFINES)
endif()

list(APPEND VIENNALS_DEPENDENCIES "viennahrle_external")

# whether to use VTK or not
option(VIENNALS_USE_VTK "Build with VTK file support." ON)
if(VIENNALS_USE_VTK)
  list(APPEND VIENNALS_DEPENDENCIES "vtk_external")
endif()

option(VIENNALS_BUILD_PYTHON_2 "Build for python2.x." OFF)
option(VIENNALS_BUILD_PYTHON_3 "Build for python3.x. Trumps python2.x build." OFF)
if(VIENNALS_BUILD_PYTHON_2 OR VIENNALS_BUILD_PYTHON_3)
  list(APPEND VIENNALS_DEPENDENCIES "pybind11_external")
endif()

option(VIENNALS_PRECOMPILE_HEADERS "Build template specialisations for shorter compile times." OFF)
# Unfortunately we cannot set this as a dependent option, since we might need to compile the
# dependencies before we can enable the PRECOMPILE_HEADERS option.
option(VIENNALS_STATIC_BUILD "Build template specialiciations in a way that allows static linking."
       OFF)

# build dependencies
set(DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)

# Include all external dependencies
include(ExternalProject)
add_custom_target(buildDependencies)
add_subdirectory(external/upstream)

# Save the dependencies in the cache (only required for python externalproject)
set(VIENNALS_DEPENDENCIES
    ${VIENNALS_DEPENDENCIES}
    CACHE INTERNAL "List of dependencies that need to be built before ViennaLS can be used" FORCE)

include(GNUInstallDirs)

set(ViennaLS_BINARY_DIR
    ${PROJECT_BINARY_DIR}
    CACHE PATH "Path to local ViennaLS installation" FORCE)

find_package(OpenMP REQUIRED)
list(APPEND VIENNALS_LIBRARIES OpenMP::OpenMP_CXX)

if(VIENNALS_USE_VTK)
  add_compile_definitions(VIENNALS_USE_VTK)

  # use VTK_USE_FILE for older versions of VTK
  if(VTK_VERSION VERSION_LESS "9.0")
    message(STATUS "Using old VTK Version: ${VTK_VERSION}")
    set(VTK_COMPONENTS
        vtkCommonCore;vtkIOCore;vtkIOXML;vtkIOXMLParser;vtkFiltersCore;vtkFiltersGeneral;vtkFiltersGeometry
    )
    # find_package is only required for populating the VTK_USE_FILE variable. This case is only
    # invoked when an external VTK installation is used since our externalproject build already
    # provides a newer version of VTK.
    find_package(VTK CONFIG REQUIRED PATHS ${VTK_DIR} COMPONENTS ${VTK_COMPONENTS} NO_DEFAULT_PATH)
    set(VTK_LIBRARIES
        vtksys;vtkIOCore;vtkexpat;vtklz4;vtkzlib;vtklzma;vtkdoubleconversion;vtkCommonMisc;vtkCommonSystem;vtkIOXML;vtkFiltersCore;vtkFiltersGeneral;vtkFiltersGeometry;vtkCommonExecutionModel;vtkCommonDataModel;vtkCommonColor;vtkIOXMLParser;vtkIOCore;vtkCommonMisc;vtkCommonComputationalGeometry;vtkCommonTransforms;vtkCommonMath
    )
    include(${VTK_USE_FILE})
  else()
    # Use the new component names instead of the old vtk* ones
    set(VTK_COMPONENTS
        CommonCore;IOCore;IOXML;IOXMLParser;FiltersCore;FiltersGeneral;FiltersGeometry)
    # For newer versions of VTK, only the target syntax is supported. This also automatically adds
    # the vtk include directories.
    set(VTK_LIBRARIES
        VTK::CommonComputationalGeometry;VTK::CommonCore;VTK::CommonDataModel;VTK::CommonExecutionModel;VTK::CommonMath;VTK::CommonMisc;VTK::CommonSystem;VTK::CommonTransforms;VTK::FiltersCore;VTK::FiltersGeneral;VTK::FiltersGeometry;VTK::IOCore;VTK::IOXML;VTK::IOXMLParser;VTK::WrappingTools
    )
  endif()
  list(APPEND VIENNALS_LIBRARIES ${VTK_LIBRARIES})
endif(VIENNALS_USE_VTK)

# set the correct paths for installation
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/")
endif()
set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}")

# This variable is only used by the example, test and precompiled library target, since those are
# compiled before the project is installed in its proper install location.
set(VIENNALS_BUILD_INCLUDE_DIRS "${${PROJECT_NAME}_SOURCE_DIR}/include")

# ##################################################################################################
# BUILD SHARED LIBRARIES
# ##################################################################################################
if(VIENNALS_PRECOMPILE_HEADERS)
  # Check if all dependencies are installed
  set(DEPENDENCIES_FOUND TRUE)
  find_package(ViennaHRLE QUIET CONFIG PATHS ${ViennaHRLE_DIR} NO_DEFAULT_PATH)
  if(NOT ViennaHRLE_FOUND)
    set(DEPENDENCIES_FOUND FALSE)
  endif()
  if(VIENNALS_USE_VTK)
    find_package(VTK QUIET PATHS ${VTK_DIR} COMPONENTS ${VTK_COMPONENTS} NO_DEFAULT_PATH)
    if(NOT VTK_FOUND)
      set(DEPENDENCIES_FOUND FALSE)
    endif()
  endif()
  # Show a message if one or more of the dependencies are missing
  if(NOT DEPENDENCIES_FOUND)
    message(
      FATAL_ERROR
        "One of the dependencies hasn't been built yet and thus can't be used by the precompiled libraries target.\n\
    Try the following steps to fix the problem:\n\
    1. Turn off precompiled libraries target by disabling the VIENNALS_PRECOMPILE_HEADERS option when configuring the project\n\
      (e.g. with the `-D VIENNALS_PRECOMPILE_HEADERS=OFF` flag)\n\
    2. Build the dependencies target: e.g. `make buildDependencies` when using Makefiles\n\
    3. Turn on precompiled libraries target again by re-enabling the VIENNALS_PRECOMPILE_HEADERS option\n\
      (e.g. with the `-D VIENNALS_PRECOMPILE_HEADERS=ON` flag)\n\
    4. Now the configure step should pick up the built dependencies and proceed as intended.")
  endif()

  # Add all specialisation source files
  file(GLOB SPECIALISATION_CPPS "lib/*.cpp")

  # Decide wheter to build the precompiled headers as a static or shared library
  if(VIENNALS_STATIC_BUILD)
    add_library(${PROJECT_NAME} STATIC ${SPECIALISATION_CPPS})

    # Find paths to all static library targets
    foreach(_lib ${VIENNALS_LIBRARIES})
      if(TARGET ${_lib})
        get_target_property(_target_type ${_lib} TYPE)
        # message(STATUS "${_lib} ${_target_type}")
        if(_target_type STREQUAL "STATIC_LIBRARY")
          list(APPEND static_lib_files $<TARGET_FILE:${_lib}>)
        elseif(_target_type STREQUAL "SHARED_LIBRARY")
          message(
            FATAL_ERROR "Not all required libraries can be statically linked!\n- Caused by ${_lib}")
        endif()
      endif()
    endforeach()

    if(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|GNU)$")
      # Create an archiver (ar) script that bundles all static library dependencies into one bundled
      # static library
      set(AR_FILE_CONTENT
          "CREATE ${CMAKE_BINARY_DIR}/libViennaLS.a\nADDLIB ${CMAKE_BINARY_DIR}/libViennaLS.a")
      foreach(_statlib ${static_lib_files})
        set(AR_FILE_CONTENT "${AR_FILE_CONTENT}\nADDLIB ${_statlib}")
      endforeach()
      set(AR_FILE_CONTENT "${AR_FILE_CONTENT}\nSAVE\nEND\n")

      file(
        GENERATE
        OUTPUT "${CMAKE_BINARY_DIR}/bundle.ar"
        CONTENT "${AR_FILE_CONTENT}")

      # Execute the archiver script after the target has been built
      add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${CMAKE_AR} -M < ${CMAKE_BINARY_DIR}/bundle.ar
        COMMENT "Bundling all static libraries into one."
        VERBATIM)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(AppleClang)$")
      find_program(libtool_bin libtool)

      # TODO: check if libtool found

      add_custom_command(
        TARGET ${PROJECT_NAME}
        POST_BUILD
        COMMAND ${libtool_bin} -static -o libViennaLS.a ${CMAKE_BINARY_DIR}/libViennaLS.a
                ${static_lib_files}
        COMMENT "Bundling all static libraries into one."
        VERBATIM)
    else()
      message(
        FATAL_ERROR
          "Building a static version of ${PROJECT_NAME} is currently not supported with the selected build tools."
      )
    endif()
  else()
    # Shared library
    add_library(${PROJECT_NAME} SHARED ${SPECIALISATION_CPPS})
  endif()

  target_link_libraries(${PROJECT_NAME} ${VIENNALS_LIBRARIES})
  # Include the project header files and ViennaHRLE
  target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include/"
                                                     ${VIENNAHRLE_INCLUDE_DIRS})
  # Also include VTK if we enabled VTK support (with VTK 9+ this is done automatically when we link
  # the vtk targets - the VTK_INCLUDE_DIRS variable is empty in that case)
  if(VIENNALS_USE_VTK)
    target_include_directories(${PROJECT_NAME} PRIVATE ${VTK_INCLUDE_DIRS})
  endif()
  set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
else()
  # By default, create header only library
  add_library(${PROJECT_NAME} INTERFACE)
endif()

# Adding the install interface generator expression makes sure that the include files are installed
# to the proper location (provided by GNUInstallDirs)
target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${VIENNALS_BUILD_INCLUDE_DIRS}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" INSTALL_DESTINATION ${LOCAL_CONFIG_DIR})

# ##################################################################################################
# BUILD EXAMPLES
# ##################################################################################################
option(VIENNALS_BUILD_EXAMPLES "Build examples." OFF)
if(VIENNALS_BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif(VIENNALS_BUILD_EXAMPLES)

# ##################################################################################################
# BUILD TESTS (Smaller examples for testing)
# ##################################################################################################
option(VIENNALS_BUILD_TESTS "Build tests." OFF)
if(VIENNALS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(Tests)
endif(VIENNALS_BUILD_TESTS)

# ##################################################################################################
# BUILD PYTHON MODULE
# ##################################################################################################
# TODO: Properly integrate the python externalproject with this cmake configuration
if(VIENNALS_BUILD_PYTHON_2 OR VIENNALS_BUILD_PYTHON_3)
  add_subdirectory(Wrapping)
endif()

# ##################################################################################################
# INSTALL
# ##################################################################################################
# install config files locations are provided by GNUInstallDirs
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE ${PROJECT_NAME}Targets.cmake
  DESTINATION ${LOCAL_CONFIG_DIR})

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})

# install include files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")
install(FILES ${HEADER_FILES} DESTINATION include)
