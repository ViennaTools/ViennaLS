project(AirGapDeposition LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaLS)

add_dependencies(ViennaLS_Examples ${PROJECT_NAME})
