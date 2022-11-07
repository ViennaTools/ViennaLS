# ##################################################################################################
# Check ViennaHRLE Dependency
# ##################################################################################################
# Check if the variable is a cache variable
get_property(
  is_cached
  CACHE ViennaHRLE_DIR
  PROPERTY TYPE)
if("${is_cached}" STREQUAL "" AND DEFINED ViennaHRLE_DIR)
  # Required for windows to convert backslashes to forward slashes
  file(TO_CMAKE_PATH "${ViennaHRLE_DIR}" ViennaHRLE_DIR)
  set(ViennaHRLE_DIR
      "${ViennaHRLE_DIR}"
      CACHE PATH "Path to ViennaHRLE installation" FORCE)
  find_package(ViennaHRLE REQUIRED PATHS ${ViennaHRLE_DIR} NO_DEFAULT_PATH)
  message(STATUS "ViennaHRLE_DIR: ${ViennaHRLE_DIR}")
  add_library(viennahrle_external INTERFACE) # dummy target
else()
  # First try to locate any local build
  find_package(ViennaHRLE QUIET PATHS "${DEPENDENCIES_DIR}/Install/viennahrle_external"
               NO_DEFAULT_PATH)
  if(NOT ViennaHRLE_FOUND)
    find_package(ViennaHRLE QUIET)
  endif()
  if(ViennaHRLE_FOUND)
    message(STATUS "Found ViennaHRLE at ${ViennaHRLE_DIR}.")
    add_library(viennahrle_external INTERFACE) # dummy target
  else()
    message(
      STATUS
        "Couldn't find ViennaHRLE. Dependency will be built with the buildDependencies target (e.g. `make buildDependencies`)."
    )
    ExternalProject_Add(
      viennahrle_external
      GIT_REPOSITORY https://github.com/ViennaTools/ViennaHRLE
      GIT_SHALLOW 1 # use --depth 1 when cloning
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DVIENNAHRLE_BUILD_EXAMPLES=OFF
                 -DVIENNAHRLE_BUILD_TESTS=OFF
      USES_TERMINAL_DOWNLOAD 1
      USES_TERMINAL_UPDATE 1
      GIT_PROGRESS 1)

    # Get install and build directory
    ExternalProject_Get_Property(viennahrle_external INSTALL_DIR BINARY_DIR)

    # Set the required variables
    set(ViennaHRLE_DIR "${INSTALL_DIR}")
    set(ViennaHRLE_INCLUDE_DIRS "${ViennaHRLE_DIR}/include")

    # add to buildDependencies target
    add_dependencies(buildDependencies viennahrle_external)

  endif()
  set(ViennaHRLE_DIR
      "${ViennaHRLE_DIR}"
      CACHE PATH "Path to ViennaHRLE installation" FORCE)
  set(ViennaHRLE_INCLUDE_DIRS
      "${ViennaHRLE_INCLUDE_DIRS}"
      CACHE PATH "Path to ViennaHRLE headers" FORCE)
endif()
