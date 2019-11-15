#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to calculate the normal vectors for
  every grid point stored in an lsDomain and outputting
  them to an explicit mesh. This also outputs the level set values
  stored in each grid point.
  \example calculateNormalVectors.cpp
*/

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double extent = 15;
  double gridDelta = 0.25;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1, lsSphere<double, D>(origin, radius))
      .apply();

  {
    lsDomain<double, D> sphere2(bounds, boundaryCons, gridDelta);
    origin[0] = -5.;
    lsMakeGeometry<double, D>(sphere2, lsSphere<double, D>(origin, radius))
        .apply();
  }

  std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1, 3).apply();

  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
            << std::endl;

  // normal vectors are only valid as long as the underlying
  // level set does not change
  lsCalculateNormalVectors<double, 3>(sphere1, true).apply();
  auto &normalVectors = sphere1.getNormalVectors();

  std::cout << "Number of Normal vectors: " << normalVectors.size()
            << std::endl;

  lsMesh mesh;
  lsToMesh<double, 3>(sphere1, mesh, true, true).apply();

  // also output LS values as scalar data
  std::vector<double> scalars;
  for (hrleConstSparseIterator<lsDomain<double, D>::DomainType> it(
           sphere1.getDomain());
       !it.isFinished(); ++it) {
    if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
      continue;

    scalars.push_back(double(it.getValue()));
  }
  mesh.insertNextScalarData(scalars, "LSValues");

  // set normal vectors as vectordata to mesh
  mesh.insertNextVectorData(normalVectors, "Normals");

  auto writer = lsVTKWriter();
  writer.setMesh(mesh);
  writer.setFileName("explicit.vtk");
  writer.apply();

  return 0;
}
