project(CGAL-Meshing LANGUAGES CXX)

add_executable(sphere "sphere.cpp")
target_link_libraries(sphere PRIVATE ViennaLS)

add_executable(trench "trench.cpp")
target_link_libraries(trench PRIVATE ViennaLS)