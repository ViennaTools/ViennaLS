cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(ViennaLS LANGUAGES CXX VERSION 2.1.0)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

option(VIENNALS_USE_VTK             "Build with VTK file support"                               ON)

option(VIENNALS_PRECOMPILE_HEADERS  "Build template specialisations for shorter compile times" OFF)
option(VIENNALS_STATIC_BUILD        "Build dependencies as static libraries"                   OFF)

option(VIENNALS_ENABLE_SANITIZER    "Enable Sanitizer for debug buidlds"                        ON)
option(VIENNALS_USE_SANITIZER       "Enable Sanitizer"                                         OFF)

option(VIENNALS_BUILD_EXAMPLES      "Build examples"                                           OFF)
option(VIENNALS_BUILD_TESTS         "Build tests"                                              OFF)

# --------------------------------------------------------------------------------------------------------
# Conditional Library options
# --------------------------------------------------------------------------------------------------------

include(CMakeDependentOption)

cmake_dependent_option(VIENNALS_BUILD_PYTHON "Build for Python (>=3)" OFF "VIENNALS_USE_VTK" OFF)

# --------------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------------

if (VIENNALS_STATIC_BUILD)
  message(STATUS "[ViennaLS] Enabling precompiled headers for static build")
  set(VIENNALS_PRECOMPILE_HEADERS ON)
  set(BUILD_SHARED_LIBS OFF)
endif()

if ((VIENNALS_ENABLE_SANITIZER AND CMAKE_BUILD_TYPE STREQUAL "Debug") OR VIENNALS_USE_SANITIZER)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory)
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=leak)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=leak)
  elseif(MSVC)
    add_compile_options(/fsanitize=address /fsanitize-address-use-after-return /fno-sanitize-address-vcasan-lib)
    add_link_options(/fsanitize=address /fsanitize-address-use-after-return /fno-sanitize-address-vcasan-lib)
  endif()

  message(STATUS "[ViennaLS] Using Sanitizer")
endif()

if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

set(CPM_USE_LOCAL_PACKAGES ON)

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)
add_library(ViennaTools::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_EXTENSIONS OFF CXX_STANDARD_REQUIRED ON)

if (MSVC)
  # https://learn.microsoft.com/cpp/c-runtime-library/math-constants
  # TODO: In case C++20 is adpoted any time soon: https://cppreference.com/w/cpp/numeric/constants

  target_compile_definitions(${PROJECT_NAME} INTERFACE _USE_MATH_DEFINES=1)
endif()

# --------------------------------------------------------------------------------------------------------
# Include directories
# --------------------------------------------------------------------------------------------------------

target_include_directories(${PROJECT_NAME} INTERFACE "include")

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("cmake/cpm.cmake")

CPMAddPackage(
  NAME            ViennaHRLE
  GIT_TAG         master
  GIT_REPOSITORY  "https://github.com/ViennaTools/ViennaHRLE"
)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX ViennaHRLE)

if (VIENNALS_USE_VTK)
  CPMAddPackage(
    NAME            VTK
    VERSION         9.3.0
    GIT_REPOSITORY  "https://gitlab.kitware.com/vtk/vtk"
    OPTIONS         "VTK_LEGACY_REMOVE ON"
                    "VTK_SMP_IMPLEMENTATION_TYPE \"OpenMP\""
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
  )

  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNALS_USE_VTK=1)
  
  target_link_libraries(${PROJECT_NAME} INTERFACE 
    VTK::CommonExecutionModel
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
    VTK::FiltersGeometry
  )

  if (NOT VTK_DIR)
    set(VTK_DIR "${VTK_BINARY_DIR}/lib")
  endif()
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Shared/Static Library
# --------------------------------------------------------------------------------------------------------

if (VIENNALS_PRECOMPILE_HEADERS)
  set(VIENNALS_LINKAGE SHARED)

  if (VIENNALS_STATIC_BUILD)
    set(VIENNALS_LINKAGE STATIC)
  endif()

  add_library(${PROJECT_NAME}Lib ${VIENNALS_LINKAGE})
  add_library(${PROJECT_NAME}::Lib ALIAS ${PROJECT_NAME}Lib)
  
  target_link_libraries(${PROJECT_NAME}Lib PUBLIC ${PROJECT_NAME})
  target_sources(${PROJECT_NAME}Lib PUBLIC "lib/specialisations.cpp")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Python Bindings
# --------------------------------------------------------------------------------------------------------

if (VIENNALS_BUILD_PYTHON)
  message(STATUS "[ViennaLS] Building Python Bindings")
  add_subdirectory(Python)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

if (VIENNALS_BUILD_EXAMPLES)
  message(STATUS "[ViennaLS] Building Examples")
  add_subdirectory(Examples)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if (VIENNALS_BUILD_TESTS)
  message(STATUS "[ViennaLS] Building Tests")

  enable_testing()
  add_subdirectory(Tests)
endif()
  
# --------------------------------------------------------------------------------------------------------
# Install Target
# --------------------------------------------------------------------------------------------------------

install(TARGETS ${PROJECT_NAME} 
        EXPORT  ${PROJECT_NAME}_Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# TOOD: This could be improved if the header files were in a sub-directory

file(GLOB headers "include/*.hpp")
install(FILES ${headers} DESTINATION include)
