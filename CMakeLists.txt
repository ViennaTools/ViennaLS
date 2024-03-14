cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
  ViennaLS
  LANGUAGES CXX
  VERSION 3.0.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNALS_USE_VTK "Build with VTK file support" ON)

option(VIENNALS_PRECOMPILE_HEADERS "Build template specialisations for shorter compile times" OFF)
option(VIENNALS_EXCLUDE_DEPENDENCIES "Exclude dependencies from install" OFF)
option(VIENNALS_STATIC_BUILD "Build dependencies as static libraries" OFF)

option(VIENNALS_ENABLE_SANITIZER "Enable Sanitizer for debug buidlds" OFF)
option(VIENNALS_USE_SANITIZER "Enable Sanitizer (Requires GCC or Clang)" OFF)

option(VIENNALS_BUILD_EXAMPLES "Build examples" OFF)
option(VIENNALS_BUILD_TESTS "Build tests" OFF)

# --------------------------------------------------------------------------------------------------------
# Conditional Library options
# --------------------------------------------------------------------------------------------------------

include(CMakeDependentOption)

cmake_dependent_option(VIENNALS_BUILD_PYTHON "Build for Python (>=3)" OFF "VIENNALS_USE_VTK" OFF)

# --------------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------------

include("cmake/openmp.cmake")

if(VIENNALS_STATIC_BUILD)
  message(STATUS "[ViennaLS] Enabling precompiled headers for static build")
  set(VIENNALS_PRECOMPILE_HEADERS ON)
  set(BUILD_SHARED_LIBS OFF)
endif()

if((VIENNALS_ENABLE_SANITIZER AND CMAKE_BUILD_TYPE STREQUAL "Debug") OR VIENNALS_USE_SANITIZER)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=memory)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=memory)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=leak)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=leak)
  endif()

  message(STATUS "[ViennaLS] Using Sanitizer")
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
  # https://learn.microsoft.com/cpp/c-runtime-library/math-constants TODO: In case C++20 is adpoted
  # any time soon: https://cppreference.com/w/cpp/numeric/constants

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

include("cmake/cpm.cmake")
include("cmake/vtk.cmake")

CPMAddPackage(
  NAME PackageProject
  VERSION 1.11.1
  GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake")

CPMFindPackage(
  NAME ViennaHRLE
  VERSION 0.4.0
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaHRLE"
  EXCLUDE_FROM_ALL ${VIENNALS_EXCLUDE_DEPENDENCIES})

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX ViennaHRLE)

if(VIENNALS_USE_VTK)
  enable_language(C)

  find_package(VTK QUIET)
  set(VIENNALS_SYSTEM_VTK ${VTK_FOUND})

  CPMFindPackage(
    NAME VTK
    # TODO: Change to 9.3.1 as soon as it's released
    # (https://github.com/ViennaTools/ViennaLS/pull/99#issuecomment-1832053715)
    GIT_TAG 99bd602bdbe8024c55e8382f7cf1013d42a14601
    VERSION 9.0.0
    GIT_REPOSITORY "https://gitlab.kitware.com/vtk/vtk"
    EXCLUDE_FROM_ALL ${VIENNALS_EXCLUDE_DEPENDENCIES}
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
# Setup Shared/Static Library
# --------------------------------------------------------------------------------------------------------

if(VIENNALS_PRECOMPILE_HEADERS)
  set(VIENNALS_LINKAGE SHARED)

  if(VIENNALS_STATIC_BUILD)
    set(VIENNALS_LINKAGE STATIC)
  endif()

  add_library(${PROJECT_NAME}_Lib ${VIENNALS_LINKAGE})
  add_library(${PROJECT_NAME}::Lib ALIAS ${PROJECT_NAME}_Lib)

  target_link_libraries(${PROJECT_NAME}_Lib PUBLIC ${PROJECT_NAME})
  target_sources(${PROJECT_NAME}_Lib PUBLIC "lib/specialisations.cpp")
  set_target_properties(${PROJECT_NAME}_Lib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
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
