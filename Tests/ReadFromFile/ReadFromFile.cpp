#include <iostream>

#include <lsFromSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to read a mesh from a file.
  \example ReadFromFile.cpp
*/

int main(int argc, char **argv) {
  if (argc == 1 || argc > 2) {
    std::cout << "Usage: <invocation> filename.vtk" << std::endl;
    return 0;
  }

  lsMesh mesh;

  lsVTKReader(mesh).readVTKLegacy(std::string(argv[1]));

  lsVTKWriter(mesh).writeVTKLegacy("test.vtk");

  return 0;
}
