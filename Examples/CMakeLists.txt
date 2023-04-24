cmake_minimum_required(VERSION 3.14)

project(ViennaLSExamples)

include("${CMAKE_SOURCE_DIR}/cmake/prepare.cmake")

# Point find_package to the binary directory instead of the install location
find_package(ViennaLS REQUIRED PATHS ${ViennaLS_BINARY_DIR} NO_DEFAULT_PATH)

list(PREPEND VIENNALS_INCLUDE_DIRS ${VIENNALS_BUILD_INCLUDE_DIRS})

# Exclude the buildExamples target from the `all` target
set_property(DIRECTORY PROPERTY EXCLUDE_FROM_ALL TRUE)

add_custom_target(buildExamples DEPENDS ${VIENNALS_DEPENDENCIES})

subdirlist(EXAMPLE_DIRS ${CMAKE_CURRENT_LIST_DIR})

# add all examples
foreach(subdir ${EXAMPLE_DIRS})
  add_subdirectory(${subdir})
endforeach()
