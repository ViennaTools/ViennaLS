cmake_minimum_required(VERSION 3.14)
project(ViennaLS_Python LANGUAGES CXX)

option(VIENNALS_USE_GPU "Enable GPU-accelerated BiCGSTAB linear solver via CUDA" OFF)

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

  if(VIENNALS_PACKAGE_PYTHON AND NOT VIENNALS_VTK_PYTHON_LIBS)
    # Do not install viennals.libs manually on macOS. This build-tree rpath
    # only lets delocate resolve VTK's @rpath dependencies during wheel repair;
    # delocate then vendors the dylibs into viennals/.dylibs and strips it.
    list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/_deps/vtk-build/lib")
  endif()
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

if(VIENNALS_VTK_RENDERING)
  vtk_module_autoinit(TARGETS ${MOD} MODULES ${VTK_LIBRARIES})
  target_compile_definitions(${MOD} PRIVATE VIENNALS_VTK_MODULE_INIT)
endif()

if(VIENNALS_USE_GPU)
  message(STATUS "[ViennaLS] Adding GPU support for 'viennals' python package")
  if(TARGET ViennaLS_GPU)
    add_dependencies(${MOD} ViennaLS_GPU)
  endif()
endif()

# Put the built .so in build/<pkg> for convenience
set_target_properties(
  ${MOD}
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
             RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}"
             ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PKG}")

# Install into site-packages/<pkg>/  (RELATIVE path!)
install(
  TARGETS ${MOD}
  LIBRARY DESTINATION ${PKG}
  RUNTIME DESTINATION ${PKG}
  ARCHIVE DESTINATION ${PKG}

# Package files
install(
  FILES ${CMAKE_SOURCE_DIR}/python/__init__.py
  DESTINATION ${PKG}

add_dependencies(${PROJECT_NAME} ${MOD})

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNALS_LIB_FOLDER ${CMAKE_BINARY_DIR}/${PKG}.libs)

if((VIENNALS_PACKAGE_PYTHON AND NOT APPLE) OR MSVC)
  # macOS wheels let delocate own dylib vendoring to avoid duplicate VTK copies.
  # Other package builds still need this staged dependency folder.
  viennacore_setup_vtk_env(${MOD} ${VIENNALS_LIB_FOLDER})

  if(VIENNALS_USE_GPU)
    foreach(_CUDA_TARGET CUDA::cudart CUDA::cusparse)
      if(TARGET ${_CUDA_TARGET})
        add_custom_command(
          TARGET ${MOD}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_directory $<TARGET_FILE_DIR:${_CUDA_TARGET}>
                  ${VIENNALS_LIB_FOLDER})
      endif()
    endforeach()
  endif()
endif()

if(NOT VIENNALS_PACKAGE_PYTHON)
  install(
    DIRECTORY ${VIENNALS_LIB_FOLDER}
    DESTINATION .
    OPTIONAL)
endif()
