project(Advection LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaLS)

add_dependencies(ViennaLS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)

add_executable(TimeIntegrationComparison TimeIntegrationComparison.cpp)
target_link_libraries(TimeIntegrationComparison PRIVATE ViennaLS)
add_dependencies(ViennaLS_Tests TimeIntegrationComparison)
add_test(NAME TimeIntegrationComparison COMMAND $<TARGET_FILE:TimeIntegrationComparison>)

add_executable(StencilLaxFriedrichsTest StencilLaxFriedrichsTest.cpp)
target_link_libraries(StencilLaxFriedrichsTest PRIVATE ViennaLS)
add_dependencies(ViennaLS_Tests StencilLaxFriedrichsTest)
add_test(NAME StencilLaxFriedrichsTest COMMAND $<TARGET_FILE:StencilLaxFriedrichsTest>)
