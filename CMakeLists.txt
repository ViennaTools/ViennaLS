cmake_minimum_required(VERSION 3.4)

project(
  "ViennaLS"
  VERSION 2.0.0)

add_definitions(-DVIENNALS_VERSION=${PROJECT_VERSION})

include(GNUInstallDirs)

# c++17 for inlince constexpr variables
SET(CMAKE_CXX_STANDARD "17")

# set default build type
SET(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# tell VS to export all symbols to its dll files
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE CACHE BOOL "Export all symbols")
endif()

# set whether to build static versions
option(VIENNALS_STATIC_BUILD "Link dependencies statically." OFF)
if(VIENNALS_STATIC_BUILD)
  if(NOT VTK_DIR)
    set(VTK_DIR $ENV{VTK_STATIC_DIR})
  endif(NOT VTK_DIR)
  list(APPEND VIENNALS_LIBRARIES "-static")
endif(VIENNALS_STATIC_BUILD)

# OpenMP Support
find_package(OpenMP REQUIRED)
list(APPEND VIENNALS_LIBRARIES OpenMP::OpenMP_CXX)

# VTK File type support
option(VIENNALS_USE_VTK "Build with VTK file support." ON)
if(VIENNALS_USE_VTK)
  if(NOT VTK_DIR)
    set(VTK_DIR $ENV{VTK_DIR})
  endif(NOT VTK_DIR)
  find_package(VTK QUIET)
  if(VTK_FOUND)
    message(STATUS "Found VTK")
    add_compile_definitions(VIENNALS_USE_VTK)
    include(${VTK_USE_FILE})
    # only link needed vtk libraries for static builds
    if(VIENNALS_STATIC_BUILD)
      set(VTK_LIBRARIES vtksys;vtkIOCore;vtkexpat;vtklz4;vtkzlib;vtklzma;vtkdoubleconversion;vtkCommonMisc;vtkCommonSystem;vtkIOXML;vtkFiltersCore;vtkFiltersGeneral;vtkFiltersGeometry;vtkCommonExecutionModel;vtkCommonDataModel;vtkCommonColor;vtkIOXMLParser;vtkIOCore;vtkCommonMisc;vtkCommonComputationalGeometry;vtkCommonTransforms;vtkCommonMath;)

      message(STATUS "Static build: Only linking ${VTK_LIBRARIES}")
      # vtksys;vtkCommonCore;vtkCommonMath;vtkCommonMisc;vtkCommonSystem;vtkCommonTransforms;vtkCommonDataModel;vtkCommonColor;vtkCommonExecutionModel;vtkCommonComputationalGeometry;vtkFiltersCore;vtkFiltersGeneral;vtkImagingCore;vtkImagingFourier;vtkFiltersStatistics;vtkFiltersExtraction;vtkInfovisCore;vtkFiltersGeometry;vtkFiltersSources;vtkRenderingCore;vtkzlib;vtkfreetype;vtkRenderingFreeType;vtkRenderingContext2D;vtkChartsCore;vtkDICOMParser;vtkdoubleconversion;vtklz4;vtklzma;vtkIOCore;vtkIOLegacy;vtkexpat;vtkIOXMLParser;vtkDomainsChemistry;vtkglew;vtkRenderingOpenGL2;vtkDomainsChemistryOpenGL2;vtkIOXML;vtkParallelCore;vtkFiltersAMR;vtkFiltersFlowPaths;vtkFiltersGeneric;vtkImagingSources;vtkFiltersHybrid;vtkFiltersHyperTree;vtkImagingGeneral;vtkFiltersImaging;vtkFiltersModeling;vtkFiltersParallel;vtkFiltersParallelImaging;vtkFiltersPoints;vtkFiltersProgrammable;vtkFiltersSMP;vtkFiltersSelection;vtkFiltersTexture;vtkFiltersTopology;verdict;vtkFiltersVerdict;vtkmetaio;vtkjpeg;vtkpng;vtktiff;vtkIOImage;vtkImagingHybrid;vtkInfovisLayout;vtkInteractionStyle;vtkImagingColor;vtkRenderingAnnotation;vtkRenderingVolume;vtkInteractionWidgets;vtkViewsCore;vtklibproj;vtkGeovisCore;vtkhdf5_src;vtkhdf5_hl_src;vtkIOAMR;vtkIOAsynchronous;vtkpugixml;vtkIOCityGML;vtkIOEnSight;vtknetcdf;vtkexodusII;vtkIOExodus;vtkgl2ps;vtkRenderingGL2PSOpenGL2;vtkIOExport;vtkIOExportOpenGL2;vtklibharu;vtkIOExportPDF;vtkIOGeometry;vtkIOImport;vtklibxml2;vtkIOInfovis;vtkIOLSDyna;vtkIOMINC;vtkogg;vtktheora;vtkIOMovie;vtkIONetCDF;vtkIOPLY;vtkjsoncpp;vtkIOParallel;vtkIOParallelXML;vtksqlite;vtkIOSQL;vtkIOSegY;vtkIOTecplotTable;vtkIOVeraOut;vtkIOVideo;vtkImagingMath;vtkImagingMorphological;vtkImagingStatistics;vtkImagingStencil;vtkInteractionImage;vtkRenderingContextOpenGL2;vtkRenderingImage;vtkRenderingLOD;vtkRenderingLabel;vtkRenderingVolumeOpenGL2;vtkViewsContext2D;vtkViewsInfovis
    endif(VIENNALS_STATIC_BUILD)
    list(APPEND VIENNALS_LIBRARIES ${VTK_LIBRARIES})
    list(APPEND VIENNALS_PYTHON_LIBRARIES ${VTK_LIBRARIES})
  else(VTK_FOUND)
    message(STATUS "No VTK install found: Building without VTK.")
    set(VIENNALS_USE_VTK OFF)
  endif(VTK_FOUND)
endif(VIENNALS_USE_VTK)


#################################################
# BUILD SHARED LIBRARIES
#################################################
option(VIENNALS_BUILD_SHARED_LIBS "Build template specialisations for shorter compile times." OFF)
if(VIENNALS_BUILD_SHARED_LIBS)
  if(NOT ViennaHRLE_DIR)
    set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
  endif(NOT ViennaHRLE_DIR)
  list(APPEND VIENNALS_LIBRARIES ViennaHRLE)
  find_package(ViennaHRLE REQUIRED)
  
  # add all specialisation source files
  file(GLOB SPECIALISATION_CPPS "lib/*.cpp")

  # build static or shared object for library
  if(VIENNALS_STATIC_BUILD)
    add_library(${PROJECT_NAME} STATIC ${SPECIALISATION_CPPS})
  else(VIENNALS_STATIC_BUILD)
    add_library(${PROJECT_NAME} SHARED ${SPECIALISATION_CPPS})
  endif(VIENNALS_STATIC_BUILD)

  target_link_libraries(${PROJECT_NAME} ${VIENNALS_LIBRARIES})
  target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include/")
  set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
else(VIENNALS_BUILD_SHARED_LIBS)
  add_library(${PROJECT_NAME} INTERFACE)
endif(VIENNALS_BUILD_SHARED_LIBS)

# set the correct paths for installation
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/")
set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}")

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
set(VIENNALS_BUILD_INCLUDE_DIRS "${${PROJECT_NAME}_SOURCE_DIR}/include")
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${VIENNALS_BUILD_INCLUDE_DIRS}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

# Run tests
#enable_testing()
#add_subdirectory(test)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                               COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${LOCAL_CONFIG_DIR})


#################################################
# BUILD EXAMPLES
#################################################
option(VIENNALS_BUILD_EXAMPLES "Build examples." OFF)
if(VIENNALS_BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif(VIENNALS_BUILD_EXAMPLES)


#################################################
# BUILD TESTS (Smaller examples for testing)
#################################################
option(VIENNALS_BUILD_TESTS "Build tests." OFF)
if(VIENNALS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(Tests)
endif(VIENNALS_BUILD_TESTS)


#################################################
# BUILD PYTHON MODULE
#################################################
add_subdirectory(Wrapping)


#################################################
# INSTALL
#################################################
# install config files
# locations are provided by GNUInstallDirs
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${LOCAL_CONFIG_DIR})


install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})

# install include files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")
install(FILES ${HEADER_FILES} DESTINATION include)
