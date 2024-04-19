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
  VERSION 2.11.1
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
# Binding macro
# --------------------------------------------------------------------------------------------------------

macro(setup_binding NAME FLAGS)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${VIENNALS_LIBRARY_OUTPUT_DIR}/${NAME})
  set(MODULE_NAME "_${NAME}")

  pybind11_add_module("${MODULE_NAME}" "pyWrap.cpp")
  add_dependencies(${PROJECT_NAME} ${MODULE_NAME})

  target_link_libraries(${MODULE_NAME} PUBLIC ViennaLS)
  target_compile_definitions(${MODULE_NAME} PRIVATE ${FLAGS} -DVIENNALS_MODULE_NAME=${MODULE_NAME})

  configure_file(__init__.py.in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init__.py)

  install(TARGETS ${MODULE_NAME} LIBRARY DESTINATION ${NAME})
  install(DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/" DESTINATION ${NAME})
endmacro()

# --------------------------------------------------------------------------------------------------------
# Setup Bindings
# --------------------------------------------------------------------------------------------------------

setup_binding(viennals2d -DVIENNALS_PYTHON_DIMENSION=2)
setup_binding(viennals3d -DVIENNALS_PYTHON_DIMENSION=3)

set(VIENNALS_LIB_FOLDER "${VIENNALS_LIBRARY_OUTPUT_DIR}/viennals.libs")

if(VIENNALS_VTK_ENV)
  setup_vtk_env(${MODULE_NAME} "${VIENNALS_LIB_FOLDER}")
  install(DIRECTORY "${VIENNALS_LIB_FOLDER}" DESTINATION .)
endif()
