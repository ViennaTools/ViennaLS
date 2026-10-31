#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

#include <vcTestAsserts.hpp>

namespace ls = viennals;

/**
  A volume mesh of two layers, material 0 from y = -20.2 nm to 3.1 nm and
  material 1 above it up to 12.3 nm, written to a VTU file and read back.
  FromVolumeMesh gives one level set per material, the first enclosing
  material 0 and the second both, whose values are the distances to the
  layers. Every second triangle runs clockwise, which the orientation of the
  elements sorts out.
*/
int main() {
  constexpr int D = 2;
  omp_set_num_threads(1);

  double gridDelta = 5e-10;
  double bounds[2 * D] = {-3.5e-8, 3.5e-8, -5e-8, 5e-8};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto domain = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Write mesh, two triangles per layer across the domain
  const double bottom = -2.02e-8, interface = 3.1e-9, top = 1.23e-8;
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (double y : {bottom, interface, top}) {
      mesh->insertNextNode({bounds[0], y, 0.});
      mesh->insertNextNode({bounds[1], y, 0.});
    }
    const std::array<unsigned, 3> triangles[4] = {
        {0, 1, 3}, {0, 2, 3}, {2, 3, 5}, {2, 4, 5}};
    for (const auto &triangle : triangles)
      mesh->insertNextTriangle(triangle);
    mesh->getCellData().insertNextScalarData({0., 0., 1., 1.}, "Material");
    ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTU, "initial.vtu").apply();
  }

  // Read mesh
  auto initialMesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::VTKReader<double>(initialMesh, "initial.vtu").apply();
  initialMesh->print();

  ls::FromVolumeMesh<double, D> reader(domain->getGrid(), initialMesh);
  reader.apply();
  auto levelSets = reader.getLevelSets();
  VC_TEST_ASSERT(levelSets.size() == 2);

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    ls::VTKWriter<double>(mesh, "LSsurface-" + std::to_string(i) + ".vtk")
        .apply();

    // the distance to the layers below the interface or the top
    const double upper = (i == 0) ? interface : top;
    double worst = 0.;
    unsigned points = 0;
    viennahrle::ConstSparseIterator<ls::Domain<double, D>::DomainType> it(
        levelSets[i]->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const double y = it.getStartIndices(1) * gridDelta;
      const double exact = std::max(y - upper, bottom - y);
      worst = std::max(worst, std::abs(it.getValue() * gridDelta - exact));
      ++points;
    }
    std::cout << "Level set " << i << ": " << points << " points, worst value "
              << worst / gridDelta << " grid spacings off" << std::endl;
    VC_TEST_ASSERT(points > 0);
    VC_TEST_ASSERT(worst < 1e-4 * gridDelta);
  }

  return 0;
}
