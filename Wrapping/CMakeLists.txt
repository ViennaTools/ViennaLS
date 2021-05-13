cmake_minimum_required(VERSION 3.4)

# lowercase viennaLS is the python export
project("ViennaLSPython"
        VERSION ${CMAKE_PROJECT_VERSION})

# include viennals
set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
set(ViennaLS_DIR "${ViennaLS_BINARY_DIR}")
set(ViennaLS_COMMON_TARGET ${PROJECT_NAME})


option(VIENNALS_BUILD_PYTHON_2 "Build for python2.x." OFF)
option(VIENNALS_BUILD_PYTHON_3 "Build for python3.x. Trumps python2.x build." OFF)

set(VIENNALS_PYTHON_SOURCE "pyWrap.cpp")
set(VIENNALS_PYTHON_MODULE_NAME "viennaLS")

message(STATUS "Pybind11: ${pybind11_DIR}")

if(VIENNALS_BUILD_PYTHON_3)
  set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version")

  find_package(ViennaHRLE REQUIRED)
  find_package(ViennaLS REQUIRED)
  find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})

  #make 3.6 python module
  # define dimension for build
  set(VIENNALS_PYTHON_MODULE_NAME_2D "${VIENNALS_PYTHON_MODULE_NAME}2d")

  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_2D} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_2D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE ${VIENNALS_PYTHON_LIBRARIES})
  #define compile dimension
  target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=2)

  # define dimension for build
  set(VIENNALS_PYTHON_MODULE_NAME_3D "${VIENNALS_PYTHON_MODULE_NAME}3d")

  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_3D} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE ${VIENNALS_PYTHON_LIBRARIES})
  #define compile dimension
  target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=3)
endif()

if(VIENNALS_BUILD_PYTHON_2 AND NOT VIENNALS_BUILD_PYTHON_3)
  set(PYBIND11_PYTHON_VERSION 2 CACHE STRING "Python version")

  # disable register warning, as Pyhton.h of Python 2 includes register variables
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-register")

  find_package(ViennaHRLE REQUIRED)
  find_package(ViennaLS REQUIRED)
  find_package(pybind11)

  #make 2.7 python module
  set(VIENNALS_PYTHON_MODULE_NAME_2D "${VIENNALS_PYTHON_MODULE_NAME}2d")

  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_2D} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_2D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE ${VIENNALS_PYTHON_LIBRARIES})
  #define compile dimension
  target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=2)

  # define dimension for build
  set(VIENNALS_PYTHON_MODULE_NAME_3D "${VIENNALS_PYTHON_MODULE_NAME}3d")

  pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_3D} ${VIENNALS_PYTHON_SOURCE})
  target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE ${VIENNALS_PYTHON_LIBRARIES})
  #define compile dimension
  target_compile_definitions(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE -DVIENNALS_PYTHON_DIMENSION=3)
endif()
