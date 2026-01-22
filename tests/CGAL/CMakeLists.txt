project(CGAL-Meshing LANGUAGES CXX)

find_package(CGAL REQUIRED)

add_executable(sphere "sphere.cpp")
target_link_libraries(sphere PRIVATE ViennaLS CGAL::CGAL)

add_executable(trench "trench.cpp")
target_link_libraries(trench PRIVATE ViennaLS CGAL::CGAL)