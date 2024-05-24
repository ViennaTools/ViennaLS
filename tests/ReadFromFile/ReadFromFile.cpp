#include <iostream>

#include <lsFromSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to read a mesh from a file.
  \example ReadFromFile.cpp
*/

namespace ls = viennals;

int main(int argc, char **argv) {
  if (argc == 1 || argc > 2) {
    std::cout << "Usage: <invocation> filename.vtk" << std::endl;
    return 0;
  }

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  ls::VTKReader<double>(mesh, std::string(argv[1])).apply();

  ls::VTKWriter<double>(mesh, "test.vtk").apply();

  return 0;
}
