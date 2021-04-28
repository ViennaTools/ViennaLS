cmake_minimum_required(VERSION 3.4)

project("VisualizationMesh")

find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

if(VIENNALS_USE_VTK)
  add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
  target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNALS_LIBRARIES})
endif()