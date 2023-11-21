cmake_minimum_required(VERSION 3.14)
project(ViennaLSPython LANGUAGES CXX)

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------

include("../cmake/cpm.cmake")

set(PYBIND11_FINDPYTHON ON)
set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version")

CPMAddPackage(
  NAME            pybind11
  VERSION         2.11.1
  GIT_REPOSITORY  "https://github.com/pybind/pybind11"
)

# --------------------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------------------

message("${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

set(VIENNALS_PYTHON_MODULE_NAME "viennals")
set(VIENNALS_LIBRARY_OUTPUT_DIR ${CMAKE_BINARY_DIR})

# --------------------------------------------------------------------------------------------------------
# Setup 2D-Bindings
# --------------------------------------------------------------------------------------------------------

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${VIENNALS_LIBRARY_OUTPUT_DIR}/viennals2d)
set(VIENNALS_PYTHON_MODULE_NAME_2D "_${VIENNALS_PYTHON_MODULE_NAME}2d")

pybind11_add_module("${VIENNALS_PYTHON_MODULE_NAME_2D}" "pyWrap.cpp")

target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_2D} PUBLIC ViennaLS)
target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=2 -DVIENNALS_MODULE_NAME=${VIENNALS_PYTHON_MODULE_NAME_2D})

set(MODULE_NAME ${VIENNALS_PYTHON_MODULE_NAME_2D})
configure_file(__init__.py.in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init__.py)

# --------------------------------------------------------------------------------------------------------
# Setup 3D-Bindings
# --------------------------------------------------------------------------------------------------------

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${VIENNALS_LIBRARY_OUTPUT_DIR}/viennals3d)
set(VIENNALS_PYTHON_MODULE_NAME_3D "_${VIENNALS_PYTHON_MODULE_NAME}3d")

pybind11_add_module("${VIENNALS_PYTHON_MODULE_NAME_3D}" "pyWrap.cpp")

target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_3D} PUBLIC ViennaLS)
target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=3 -DVIENNALS_MODULE_NAME=${VIENNALS_PYTHON_MODULE_NAME_3D})

set(MODULE_NAME ${VIENNALS_PYTHON_MODULE_NAME_3D})
configure_file(__init__.py.in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init__.py)
