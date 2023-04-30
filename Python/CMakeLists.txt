cmake_minimum_required(VERSION 3.14)

project(ViennaLSPython)

include("${CMAKE_SOURCE_DIR}/cmake/prepare.cmake")

# Point find_package to the binary directory instead of the install location
find_package(ViennaLS REQUIRED PATHS ${ViennaLS_BINARY_DIR} NO_DEFAULT_PATH)

list(PREPEND VIENNALS_INCLUDE_DIRS ${VIENNALS_BUILD_INCLUDE_DIRS})

message(STATUS "ViennaLS version: ${VIENNALS_VERSION}")

set(PYBIND11_PYTHON_VERSION
    3
    CACHE STRING "Python version")

find_package(pybind11 REQUIRED PATHS ${pybind11_DIR} NO_DEFAULT_PATH)

set(VIENNALS_PYTHON_SOURCE pyWrap.cpp)
set(VIENNALS_PYTHON_MODULE_NAME "viennals")

# ##################################################################################################
# BUILD 2D PYTHON LIBRARY
# ##################################################################################################
set(VIENNALS_PYTHON_MODULE_NAME_2D "_${VIENNALS_PYTHON_MODULE_NAME}2d")

pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_2D} ${VIENNALS_PYTHON_SOURCE})
target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_2D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_2D} PRIVATE ${VIENNALS_LIBRARIES})
target_compile_definitions(
  ${VIENNALS_PYTHON_MODULE_NAME_2D}
  PRIVATE -DVIENNALS_PYTHON_DIMENSION=2 -DVIENNALS_MODULE_NAME=${VIENNALS_PYTHON_MODULE_NAME_2D})

message(STATUS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY_Release})
message(STATUS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY_Debug})
message(STATUS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
message(STATUS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_Release})
message(STATUS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_Debug})

if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(MODULE_NAME ${VIENNALS_PYTHON_MODULE_NAME_2D})
  if(WIN32)
    set_target_properties(
      ${VIENNALS_PYTHON_MODULE_NAME_2D} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                   ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals2d)
  else()
    set_target_properties(
      ${VIENNALS_PYTHON_MODULE_NAME_2D} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                   ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals2d)
  endif()

  configure_file(__init__.py.in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals2d/__init__.py)
endif()

# ##################################################################################################
# BUILD 3D PYTHON LIBRARY
# ##################################################################################################
set(VIENNALS_PYTHON_MODULE_NAME_3D "_${VIENNALS_PYTHON_MODULE_NAME}3d")

pybind11_add_module(${VIENNALS_PYTHON_MODULE_NAME_3D} ${VIENNALS_PYTHON_SOURCE})
target_include_directories(${VIENNALS_PYTHON_MODULE_NAME_3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${VIENNALS_PYTHON_MODULE_NAME_3D} PRIVATE ${VIENNALS_LIBRARIES})
target_compile_definitions(
  ${VIENNALS_PYTHON_MODULE_NAME_3D}
  PRIVATE -DVIENNALS_PYTHON_DIMENSION=3 -DVIENNALS_MODULE_NAME=${VIENNALS_PYTHON_MODULE_NAME_3D})

if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(MODULE_NAME ${VIENNALS_PYTHON_MODULE_NAME_3D})

  if(WIN32)
    set_target_properties(
      ${VIENNALS_PYTHON_MODULE_NAME_3D} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                   ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals3d)
  else()
    set_target_properties(
      ${VIENNALS_PYTHON_MODULE_NAME_3D} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                                   ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals3d)
  endif()

  configure_file(__init__.py.in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/viennals3d/__init__.py)
endif()
