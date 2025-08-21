cmake_minimum_required(VERSION 3.14)
project(ViennaLS_Python LANGUAGES CXX)

add_custom_target(${PROJECT_NAME} ALL)

# --------------------------------------------------------------------------------------------------------
# Global CMake Configuration
# └ See: https://github.com/ViennaTools/ViennaPS/blob/c76e371817a797dfe2800691f00cb93317b731fa/python/CMakeLists.txt#L8
# --------------------------------------------------------------------------------------------------------

set(CMAKE_MACOSX_RPATH ON)

set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

if(NOT APPLE)
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN")
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../vtkmodules")
  list(APPEND CMAKE_INSTALL_RPATH "$ORIGIN/../viennals.libs")
else()
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path")
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../vtkmodules")
  list(APPEND CMAKE_INSTALL_RPATH "@loader_path/../viennals.libs")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("../cmake/cpm.cmake")

set(PYBIND11_FINDPYTHON ON)

CPMFindPackage(
  NAME pybind11
  VERSION 3.0.0
  GIT_REPOSITORY "https://github.com/pybind/pybind11")

# --------------------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------------------

set(VIENNALS_PYTHON_MODULE_NAME "viennals")
set(VIENNALS_LIBRARY_OUTPUT_DIR ${CMAKE_BINARY_DIR})

if(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(VIENNALS_LIBRARY_OUTPUT_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
  message(STATUS "[ViennaLS] Using Library Output Directory for bindings")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Bindings
# --------------------------------------------------------------------------------------------------------

viennacore_setup_binding(viennals ${VIENNALS_LIBRARY_OUTPUT_DIR})

target_link_libraries(viennals PRIVATE ViennaLS)

target_compile_definitions(
  viennals PRIVATE -DVIENNALS_MODULE_NAME=viennals
                   -DVIENNALS_VERSION=${CMAKE_PROJECT_VERSION})

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNALS_LIB_FOLDER ${VIENNALS_LIBRARY_OUTPUT_DIR}/viennals.libs)

# Not required for both targets, one will suffice
viennacore_setup_vtk_env(viennals ${VIENNALS_LIB_FOLDER})

install(
  DIRECTORY ${VIENNALS_LIB_FOLDER}
  DESTINATION .
  OPTIONAL)
