#include <cmath>
#include <iostream>
#include <string>

#include <lsMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

#include <vcTestAsserts.hpp>

/**
  A mesh of doubles written by VTKWriter and read back by VTKReader keeps its
  coordinates: a Euclidean surface is exact to rounding, and single precision
  would move a point at 15.9 by up to 5e-7. Writes and reads a triangle with
  such coordinates as polydata (.vtp) and as unstructured grid (.vtu).
  \example VTKPrecision.cpp
*/

namespace ls = viennals;

int main() {
  const double coordinates[3][3] = {{15.874386314, 1. / 3., -2. / 7.},
                                    {-7.123456789012, 0.1, 3.3},
                                    {1e-9, 22.000000123, -0.7}};
  for (const std::string extension : {".vtp", ".vtu"}) {
    auto mesh = ls::SmartPointer<ls::Mesh<double>>::New();
    for (const auto &c : coordinates)
      mesh->insertNextNode({c[0], c[1], c[2]});
    mesh->insertNextTriangle({0, 1, 2});
    const std::string file = "VTKPrecision" + extension;
    ls::VTKWriter<double>(mesh, file).apply();
    auto read = ls::SmartPointer<ls::Mesh<double>>::New();
    ls::VTKReader<double>(read, file).apply();
    VC_TEST_ASSERT(read->getNodes().size() == 3);
    double largest = 0.;
    for (unsigned n = 0; n < 3; ++n)
      for (int i = 0; i < 3; ++i)
        largest = std::max(
            largest, std::abs(read->getNodes()[n][i] - coordinates[n][i]));
    std::cout << extension << ": largest change of a coordinate " << largest
              << std::endl;
    VC_TEST_ASSERT(largest < 1e-12);
  }
  return 0;
}
