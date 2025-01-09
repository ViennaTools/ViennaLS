#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMeshRefined.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write different
  meshes created by the algorithms lsToVoxelMesh and lsToSurfaceMesh.
  \example Make3DSphere.cpp
*/

namespace ls = viennals;

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double gridDelta = 0.4;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);

  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);
  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = -5.;
  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  ls::BooleanOperation<double, D>(sphere1, sphere2,
                                  ls::BooleanOperationEnum::UNION)
      .apply();

  std::cout << "Number of points: " << sphere1->getDomain().getNumberOfPoints()
            << std::endl;
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  ls::ToSurfaceMeshRefined<double, double, D>(sphere1, mesh).apply();
  ls::VTKWriter<double>(mesh, "test-refined.vtp").apply();

  std::cout << "Refined mesh written to test-refined.vtp" << std::endl;
  std::cout << "Number of points: " << mesh->nodes.size() << std::endl;
  std::cout << "Number of triangles: " << mesh->triangles.size() << std::endl;

  auto &nodes = mesh->nodes;
  auto minNodeDistance =
      ls::ToSurfaceMeshRefined<double, double, D>::minNodeDistance;
  for (unsigned i = 0; i < nodes.size(); ++i) {
    for (unsigned j = i + 1; j < nodes.size(); ++j) {
      double dist = 0.;
      // manhattan distance
      for (unsigned k = 0; k < D; ++k) {
        dist += std::abs(nodes[i][k] - nodes[j][k]);
      }
      if (dist < minNodeDistance) {
        std::cout << "Distance between nodes " << i << " and " << j
                  << " is smaller than minNodeDistance: " << dist << std::endl;
      }
    }
  }
  std::cout << "Minimum node distance: " << minNodeDistance << std::endl;

  ls::ToSurfaceMeshRefined<double, double, D> noCheck(sphere1, mesh);
  noCheck.setCheckNodeDistance(false);
  noCheck.apply();

  std::cout << "Not refined mesh written to test-not-refined.vtp" << std::endl;
  std::cout << "Number of points: " << mesh->nodes.size() << std::endl;
  std::cout << "Number of triangles: " << mesh->triangles.size() << std::endl;

  ls::VTKWriter<double>(mesh, "test-no-refined.vtp").apply();

  return 0;
}
