project(Lithography LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaLS)

configure_file(polygon.csv ${CMAKE_CURRENT_BINARY_DIR}/polygon.csv COPYONLY)
