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
# Setup Bindings
# --------------------------------------------------------------------------------------------------------

set(PKG viennals)
set(MOD _core)

pybind11_add_module(${MOD} pyWrap.cpp)
target_link_libraries(${MOD} PRIVATE ViennaLS)
target_compile_definitions(${MOD} PRIVATE VIENNALS_MODULE_NAME=${MOD})
target_compile_definitions(${MOD} PRIVATE VIENNALS_VERSION=${PROJECT_VERSION})

# Put the built .so in build/<pkg> for convenience
set_target_properties(${MOD} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
)

# Install into site-packages/<pkg>/  (RELATIVE path!)
install(TARGETS ${MOD}
  LIBRARY DESTINATION ${PKG}
  RUNTIME DESTINATION ${PKG}
  ARCHIVE DESTINATION ${PKG}
)

# Package files 
install(FILES ${CMAKE_SOURCE_DIR}/python/__init__.py DESTINATION ${PKG})

add_dependencies(${PROJECT_NAME} ${MOD})

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNALS_LIB_FOLDER ${CMAKE_BINARY_DIR}/${PKG}.libs)

# Not required for both targets, one will suffice
viennacore_setup_vtk_env(${MOD} ${VIENNALS_LIB_FOLDER})

install(
  DIRECTORY ${VIENNALS_LIB_FOLDER}
  DESTINATION .
  OPTIONAL)
