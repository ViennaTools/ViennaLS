cmake_minimum_required(VERSION 3.8)

project("SmartPointer")

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNALS_LIBRARIES})

add_dependencies(buildTests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
