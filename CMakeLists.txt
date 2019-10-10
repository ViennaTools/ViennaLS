cmake_minimum_required(VERSION 3.4)

project(
  "ViennaLS"
  VERSION 0.0.1)

include(GNUInstallDirs)

# needed because g++ optimizes away the T(constexpr)
# workaround for passing references
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_CXX_STANDARD "17")
endif()

#################################################
# BUILD SHARED LIBRARIES
#################################################
option(VIENNALS_BUILD_SHARED_LIBS "Build template specialisations for shorter compile times." ON)
if(VIENNALS_BUILD_SHARED_LIBS)
  # add all specialisation source files
  set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
  find_package(ViennaHRLE REQUIRED)
  file(GLOB SPECIALISATION_CPPS "specialisations/*.cpp")
  add_library(${PROJECT_NAME} SHARED ${SPECIALISATION_CPPS})
  target_link_libraries(${PROJECT_NAME} ViennaHRLE)
  target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include/")
  set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
else(VIENNALS_BUILD_SHARED_LIBS)
  add_library(${PROJECT_NAME} INTERFACE)
endif(VIENNALS_BUILD_SHARED_LIBS)

# set the correct paths for installation
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}-${PROJECT_VERSION}/")
set(LOCAL_CONFIG_DIR "lib/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
set(VIENNALS_BUILD_INCLUDE_DIRS "${${PROJECT_NAME}_SOURCE_DIR}/include" "${${PROJECT_NAME}_SOURCE_DIR}/specialisations")
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${VIENNALS_BUILD_INCLUDE_DIRS}>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

# Run tests
#enable_testing()
#add_subdirectory(test)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                               COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${LOCAL_CONFIG_DIR})


#################################################
# BUILD EXAMPLES
#################################################
option(VIENNALS_BUILD_EXAMPLES "Build examples." OFF)
if(VIENNALS_BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif(VIENNALS_BUILD_EXAMPLES)


#################################################
# INSTALL
#################################################
# install config files
# locations are provided by GNUInstallDirs
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(EXPORT ${PROJECT_NAME}_Targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${LOCAL_CONFIG_DIR})


install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${LOCAL_CONFIG_DIR})

# install include files
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.hpp")
install(FILES ${HEADER_FILES} DESTINATION include)
file(GLOB_RECURSE HEADER_FILES "${PROJECT_SOURCE_DIR}/specialisations/*.hpp")
install(FILES ${HEADER_FILES} DESTINATION include)
