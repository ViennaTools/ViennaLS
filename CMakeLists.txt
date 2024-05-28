cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
  ViennaLS
  LANGUAGES CXX
  VERSION 3.2.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNALS_USE_VTK "Build with VTK file support" ON)

option(VIENNALS_PRECOMPILE_HEADERS "Build template specialisations for shorter compile times" OFF)
option(VIENNALS_STATIC_BUILD "Build dependencies as static libraries" OFF)
option(VIENNALS_ENABLE_SANITIZER "Enable Sanitizers" OFF)

option(VIENNALS_VTK_PYTHON_LIBS "Build against shared libraries shipped by VTKs Python library" OFF)
option(VIENNALS_PACKAGE_PYTHON "Build python bindings with intent to publish wheel" OFF)
option(VIENNALS_BUILD_PYTHON "Build python bindings" OFF)

option(VIENNALS_BUILD_EXAMPLES "Build examples" OFF)
option(VIENNALS_BUILD_TESTS "Build tests" OFF)

# --------------------------------------------------------------------------------------------------------
# Global CMake Configuration
# └ See: https://github.com/ViennaTools/ViennaPS/blob/c76e371817a797dfe2800691f00cb93317b731fa/CMakeLists.txt#L30
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_VTK_PYTHON_LIBS AND VIENNALS_PACKAGE_PYTHON)
  message(
    STATUS "[ViennaLS] Disabling `VIENNALS_PACKAGE_PYTHON` in favor of `VIENNALS_VTK_PYTHON_LIBS`")
  set(VIENNALS_PACKAGE_PYTHON OFF)
endif()

if(VIENNALS_PACKAGE_PYTHON)
  message(
    STATUS
      "[ViennaLS] Using remote VTK due to 'VIENNALS_PACKAGE_PYTHON', the build will take a long time!"
  )

  # CIBUILDWHEEL is absolutely cursed and won't work when this is turned on.
  # I don't know why...
  if(NOT VIENNALS_IS_CI)
    set(CPM_DOWNLOAD_ALL ON)
  endif()

  set(VIENNALS_BUILD_PYTHON ON)
endif()

if(VIENNALS_BUILD_PYTHON)
  set(CMAKE_PLATFORM_NO_VERSIONED_SONAME ON)
endif()

# --------------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_BUILD_PYTHON)
  message(STATUS "[ViennaLS] Using VTK for python bindings")
  set(VIENNALS_USE_VTK ON)
endif()

if(VIENNALS_STATIC_BUILD)
  message(STATUS "[ViennaLS] Enabling precompiled headers for static build")
  set(VIENNALS_PRECOMPILE_HEADERS ON)
  set(BUILD_SHARED_LIBS OFF)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)
add_library(ViennaTools::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 17
             CXX_EXTENSIONS OFF
             CXX_STANDARD_REQUIRED ON)

if(MSVC)
  # https://learn.microsoft.com/cpp/c-runtime-library/math-constants
  # TODO: In case C++20 is adpoted any time soon: https://cppreference.com/w/cpp/numeric/constants
  target_compile_definitions(${PROJECT_NAME} INTERFACE _USE_MATH_DEFINES=1)
endif()

# --------------------------------------------------------------------------------------------------------
# Include directories
# --------------------------------------------------------------------------------------------------------

target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/viennals>
                            $<INSTALL_INTERFACE:include/viennals-${PROJECT_VERSION}>)

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include(cmake/cpm.cmake)
include(cmake/vtk.cmake)

CPMAddPackage(
  NAME ViennaCore
  GIT_TAG main # TODO: Create Tag
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaCore"
  OPTIONS "VIENNACORE_FORMAT_EXCLUDE docs/")

CPMAddPackage(
  NAME PackageProject
  VERSION 1.11.1
  GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake"
  EXCLUDE_FROM_ALL ${VIENNALS_BUILD_PYTHON})

CPMFindPackage(
  NAME ViennaHRLE
  VERSION 0.4.0
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaHRLE"
  EXCLUDE_FROM_ALL ${VIENNALS_BUILD_PYTHON})

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX ViennaHRLE ViennaCore)

if(VIENNALS_USE_VTK AND VIENNALS_VTK_PYTHON_LIBS)
  import_vtk_python()
  target_link_libraries(${PROJECT_NAME} INTERFACE VTK::PythonLibs)
elseif(VIENNALS_USE_VTK)
  enable_language(C)

  find_package(VTK QUIET)
  set(VIENNALS_SYSTEM_VTK
      ${VTK_FOUND}
      CACHE INTERNAL "")

  CPMFindPackage(
    NAME VTK
    # TODO: Change to 9.3.1 as soon as it's released
    # (https://github.com/ViennaTools/ViennaLS/pull/99#issuecomment-1832053715)
    GIT_TAG 99bd602bdbe8024c55e8382f7cf1013d42a14601
    VERSION 9.0.0
    GIT_REPOSITORY "https://gitlab.kitware.com/vtk/vtk"
    EXCLUDE_FROM_ALL ${VIENNALS_BUILD_PYTHON}
    OPTIONS "VTK_LEGACY_REMOVE ON"
            "VTK_SMP_IMPLEMENTATION_TYPE \"OpenMP\""
            "VTK_GROUP_ENABLE_Rendering NO"
            "VTK_GROUP_ENABLE_Imaging NO"
            "VTK_GROUP_ENABLE_Views NO"
            "TK_GROUP_ENABLE_Web NO"
            "VTK_GROUP_ENABLE_Qt NO"
            "VTK_GROUP_ENABLE_MPI NO"
            "VTK_ENABLE_WRAPPING NO"
            "VTK_MODULE_ENABLE_VTK_libproj NO"
            "VTK_MODULE_ENABLE_VTK_CommonExecutionModel YES"
            "VTK_MODULE_ENABLE_VTK_CommonMisc YES"
            "VTK_MODULE_ENABLE_VTK_CommonSystem YES"
            "VTK_MODULE_ENABLE_VTK_CommonMath YES"
            "VTK_MODULE_ENABLE_VTK_CommonCore YES"
            "VTK_MODULE_ENABLE_VTK_CommonTransforms YES"
            "VTK_MODULE_ENABLE_VTK_CommonComputationalGeometry YES"
            "VTK_MODULE_ENABLE_VTK_IOCore YES"
            "VTK_MODULE_ENABLE_VTK_IOXMLParser YES"
            "VTK_MODULE_ENABLE_VTK_IOXML YES"
            "VTK_MODULE_ENABLE_VTK_FiltersCore YES"
            "VTK_MODULE_ENABLE_VTK_FiltersGeneral YES"
            "VTK_MODULE_ENABLE_VTK_FiltersGeometry YES"
            # This fixes an issue with the install script that causes find_package to fail.
            "CMAKE_PROJECT_NAME VTK")

  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNALS_USE_VTK=1)

  target_link_libraries(
    ${PROJECT_NAME}
    INTERFACE VTK::CommonExecutionModel
              VTK::CommonMisc
              VTK::CommonSystem
              VTK::CommonMath
              VTK::CommonCore
              VTK::CommonTransforms
              VTK::CommonComputationalGeometry
              VTK::IOCore
              VTK::IOXMLParser
              VTK::IOXML
              VTK::FiltersCore
              VTK::FiltersGeneral
              VTK::FiltersGeometry)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Sanitizer
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_ENABLE_SANITIZER)
  viennacore_enable_sanitizer()
endif()

# --------------------------------------------------------------------------------------------------------
# Check need for environment setups
# --------------------------------------------------------------------------------------------------------

if(NOT VIENNALS_VTK_PYTHON_LIBS AND (NOT VIENNALS_SYSTEM_VTK OR VIENNALS_PACKAGE_PYTHON))
  set(VIENNALS_VTK_ENV
      ON
      CACHE INTERNAL "")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Shared/Static Library
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_PRECOMPILE_HEADERS)
  set(VIENNALS_LINKAGE SHARED)

  if(VIENNALS_STATIC_BUILD)
    set(VIENNALS_LINKAGE STATIC)
  endif()

  add_library(${PROJECT_NAME}Lib ${VIENNALS_LINKAGE})
  add_library(${PROJECT_NAME}::Lib ALIAS ${PROJECT_NAME}Lib)

  target_link_libraries(${PROJECT_NAME}Lib PUBLIC ${PROJECT_NAME})
  target_sources(${PROJECT_NAME}Lib PUBLIC "lib/specialisations.cpp")
  set_target_properties(${PROJECT_NAME}Lib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_BUILD_EXAMPLES)
  message(STATUS "[ViennaLS] Building Examples")
  add_subdirectory(examples)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_BUILD_TESTS)
  message(STATUS "[ViennaLS] Building Tests")

  enable_testing()
  add_subdirectory(tests)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Python Bindings
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_BUILD_PYTHON)
  message(STATUS "[ViennaLS] Building Python Bindings")
  add_subdirectory(python)
  return()
endif()

# --------------------------------------------------------------------------------------------------------
# Install Target
# --------------------------------------------------------------------------------------------------------

packageProject(
  NAME ${PROJECT_NAME} NAMESPACE ViennaTools
  VERSION ${PROJECT_VERSION}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/viennals
  INCLUDE_DESTINATION include/viennals-${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  DEPENDENCIES "VTK 9.0.0;ViennaHRLE 0.4.0")
