#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double extent = 15;
  double gridDelta = 0.25;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain_double_3::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain_double_3::BoundaryType::SYMMETRIC_BOUNDARY;
  lsDomain_double_3 sphere1(bounds, boundaryCons, gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius);

  {
    lsDomain_double_3 sphere2(bounds, boundaryCons, gridDelta);
    origin[0] = -5.;
    lsMakeGeometry<double, D>(sphere2).makeSphere(origin, radius);
  }

  std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1).apply(3);

  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
            << std::endl;

  std::vector<hrleVectorType<double, 3>> normalVectors;
  lsCalculateNormalVectors<double, 3>(sphere1, normalVectors).apply();

  std::cout << "Number of Normal vectors: " << normalVectors.size()
            << std::endl;

  lsMesh mesh;
  lsToMesh<double, 3>(sphere1, mesh, true).apply();

  // also output LS values as scalar data
  std::vector<double> scalars;
  for (hrleConstSparseIterator<lsDomain_double_3::DomainType> it(
           sphere1.getDomain());
       !it.isFinished(); ++it) {
    if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
      continue;

    scalars.push_back(double(it.getValue()));
  }
  mesh.insertNextScalarData(scalars, "LSValues");

  // set normal vectors as vectordata to mesh
  mesh.insertNextVectorData(normalVectors, "Normals");

  auto writer = lsVTKWriter(mesh);
  writer.writeVTKLegacy("explicit.vtk");

  return 0;
}
