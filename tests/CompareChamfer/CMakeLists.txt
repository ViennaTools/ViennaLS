project(CompareChamfer LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaLS)

# Add OpenMP support for MSVC
if(MSVC)
  target_compile_options(${PROJECT_NAME} PRIVATE /openmp)
endif()

add_dependencies(ViennaLS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
