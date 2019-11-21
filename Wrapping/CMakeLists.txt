cmake_minimum_required(VERSION 3.4)

# lowercase viennaLS is the python export
project(ViennaLSPython)

# include viennals
set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
set(ViennaLS_DIR "${ViennaLS_BINARY_DIR}")
set(ViennaLS_COMMON_TARGET ${PROJECT_NAME})
message(STATUS "Dir: ${ViennaLS_BINARY_DIR}")
find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

option(VIENNALS_BUILD_PYTHON_2_7 "Build for python2.7." OFF)
option(VIENNALS_BUILD_PYTHON_3_6 "Build for python3.6. Trumps python2.7 build." OFF)

set(VIENNALS_PYTHON_SOURCE "pyWrap.cpp")
set(VIENNALS_PYTHON_MODULE_NAME "viennaLS")

# set the correct path for pybind11
if(NOT pybind11_DIR)
  set(pybind11_DIR $ENV{PYBIND11_DIR})
endif(NOT pybind11_DIR)

if(VIENNALS_BUILD_PYTHON_3_6)
  set(PYBIND11_PYTHON_VERSION 3.6 CACHE STRING "Python version")

  find_package(pybind11)

  #make 2.7 python module
  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME} PRIVATE ${VIENNALS_LIBRARIES})
endif()

if(VIENNALS_BUILD_PYTHON_2_7 AND NOT EXAMPLE_BUILD_PYTHON_3_6)
  set(PYBIND11_PYTHON_VERSION 2.7 CACHE STRING "Python version")

  find_package(pybind11)

  #make 2.7 python module
  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME} PRIVATE ${VIENNALS_LIBRARIES})
endif()
