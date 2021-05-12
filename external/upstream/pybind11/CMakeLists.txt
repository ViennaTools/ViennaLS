externalproject_add(pybind11-external
  GIT_REPOSITORY
  https://github.com/pybind/pybind11.git
  # Hash of tags/v2.6.2
  GIT_TAG #v2.6.2
  8de7772cc72daca8e947b79b83fea46214931604
  # GIT_SHALLOW 1
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    # -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    #-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    # -DBUILD_TESTING=OFF
    # -DBUILD_EXAMPLES=OFF
  # BUILD_COMMAND
  "cmake --build build -- -j4"
  # # No install needed
  # INSTALL_COMMAND
  "cmake --install build"
  UPDATE_COMMAND
    ""
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
)

# Get install and build directory
externalproject_get_property(pybind11-external INSTALL_DIR BINARY_DIR)

set(
  pybind11_DIR "${INSTALL_DIR}"
  CACHE PATH "Path to pybind11 installation"
  )

