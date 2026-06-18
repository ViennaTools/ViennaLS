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
  ARCHIVE DESTINATION ${PKG})

# Package files
install(
  FILES ${CMAKE_SOURCE_DIR}/python/__init__.py
  DESTINATION ${PKG})

add_dependencies(${PROJECT_NAME} ${MOD})

# --------------------------------------------------------------------------------------------------------
# Setup Lib-Folder
# --------------------------------------------------------------------------------------------------------

set(VIENNALS_LIB_FOLDER ${CMAKE_BINARY_DIR}/${PKG}.libs)

if((VIENNALS_PACKAGE_PYTHON AND NOT APPLE) OR MSVC)
  # macOS wheels let delocate own dylib vendoring to avoid duplicate VTK copies.
  # Other package builds still need this staged dependency folder.
  viennacore_setup_vtk_env(${MOD} ${VIENNALS_LIB_FOLDER})

  if(MSVC)
    # /openmp:llvm links against the LLVM OpenMP runtime, which is usually
    # available on GitHub runners but not on end-user Windows machines.
    get_filename_component(_VIENNALS_CLANG_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    set(_VIENNALS_OPENMP_DLL_CANDIDATES)

    if(DEFINED ENV{VCINSTALLDIR})
      list(
        APPEND
        _VIENNALS_OPENMP_DLL_CANDIDATES
        "$ENV{VCINSTALLDIR}/Tools/Llvm/bin/libomp140.x86_64.dll"
        "$ENV{VCINSTALLDIR}/Tools/Llvm/x64/bin/libomp140.x86_64.dll")
    endif()

    if(DEFINED ENV{VCToolsInstallDir})
      list(
        APPEND
        _VIENNALS_OPENMP_DLL_CANDIDATES
        "$ENV{VCToolsInstallDir}/../../Llvm/bin/libomp140.x86_64.dll"
        "$ENV{VCToolsInstallDir}/../../Llvm/x64/bin/libomp140.x86_64.dll")
    endif()

    list(
      APPEND
      _VIENNALS_OPENMP_DLL_CANDIDATES
      "${_VIENNALS_CLANG_BIN_DIR}/../../../../../../Tools/Llvm/bin/libomp140.x86_64.dll"
      "${_VIENNALS_CLANG_BIN_DIR}/../../../../../../Tools/Llvm/x64/bin/libomp140.x86_64.dll")

    unset(_VIENNALS_OPENMP_DLL)
    foreach(_VIENNALS_OPENMP_DLL_CANDIDATE IN LISTS _VIENNALS_OPENMP_DLL_CANDIDATES)
      get_filename_component(
        _VIENNALS_OPENMP_DLL_CANDIDATE_REAL "${_VIENNALS_OPENMP_DLL_CANDIDATE}" ABSOLUTE)
      if(EXISTS "${_VIENNALS_OPENMP_DLL_CANDIDATE_REAL}")
        set(_VIENNALS_OPENMP_DLL "${_VIENNALS_OPENMP_DLL_CANDIDATE_REAL}")
        break()
      endif()
    endforeach()

    if(_VIENNALS_OPENMP_DLL)
      add_custom_command(
        TARGET ${MOD}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${VIENNALS_LIB_FOLDER}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_VIENNALS_OPENMP_DLL}"
                "${VIENNALS_LIB_FOLDER}/")
    else()
      message(
        WARNING
          "[ViennaLS] Could not locate libomp140.x86_64.dll. Windows wheels built "
          "with /openmp:llvm may fail to import on machines without the LLVM OpenMP runtime."
      )
    endif()
  endif()

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
  # Dev/test builds: install the whole lib folder (VTK DLLs from vcpkg on
  # Windows, VTK .so on Linux, etc.) so the module can be imported in-place.
  install(
    DIRECTORY ${VIENNALS_LIB_FOLDER}
    DESTINATION .
    OPTIONAL)
elseif(MSVC)
  # Package build on Windows: static VTK is already compiled into _core.pyd.
  # Only DLLs are needed at runtime (libomp140.x86_64.dll etc.).
  # Installing the full VIENNALS_LIB_FOLDER would add hundreds of MBs of .lib
  # files; restrict to *.dll only.
  install(
    DIRECTORY ${VIENNALS_LIB_FOLDER}
    DESTINATION .
    OPTIONAL
    FILES_MATCHING PATTERN "*.dll")
endif()
