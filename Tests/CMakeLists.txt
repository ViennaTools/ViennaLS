cmake_minimum_required(VERSION 3.8)

project(ViennaLSTests)

include("${CMAKE_SOURCE_DIR}/cmake/prepare.cmake")

# Point find_package to the binary directory instead of the install location
find_package(ViennaLS REQUIRED CONFIG PATHS ${ViennaLS_BINARY_DIR} NO_DEFAULT_PATH)

list(PREPEND VIENNALS_INCLUDE_DIRS ${VIENNALS_BUILD_INCLUDE_DIRS})

add_custom_target(buildTests DEPENDS ${VIENNALS_DEPENDENCIES})

subdirlist(TEST_DIRS ${CMAKE_CURRENT_LIST_DIR})

# add all tests
foreach(subdir ${TEST_DIRS})
  add_subdirectory(${subdir})
endforeach()
