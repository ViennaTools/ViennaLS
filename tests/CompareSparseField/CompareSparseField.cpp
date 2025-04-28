#include <chrono>
#include <cmath>
#include <iostream>

#include <lsCompareNarrowBand.hpp>
#include <lsCompareSparseField.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsReduce.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsCompareSparseField that compares two level sets by measuring
  the SDF value differences using sparse iteration. This test creates two
  spheres with different centers and measures their differences efficiently.
  \example CompareSparseField.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  // Create first sphere (target)
  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin1[D] = {0., 0.};
  double radius1 = 5.0;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  // Create second sphere (sample) with different center
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin2[D] = {2., 1.}; // Shifted center
  // double origin2[D] = {0., 0.}; // Same center
  double radius2 = 5.0; // Same radius

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Make sure target is expanded
  ls::Expand<double, D>(sphere1, 50).apply();

  // Reduce the sample level set to a sparse field
  ls::Reduce<double, D>(sphere2, 1).apply();

  // Export both spheres as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere1_expanded_target.vtp").apply();
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere1, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere1_surface.vtp").apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere2, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere2_sparse_sample.vtp").apply();
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere2, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere2_surface.vtp").apply();
  }

  // Compare using sparse field comparison
  ls::CompareSparseField<double, D> compareSparseField(sphere1, sphere2);

  // Create mesh for visualization of differences
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  compareSparseField.setFillSampleWithDistances(true);
  compareSparseField.setOutputMesh(mesh);
  compareSparseField.apply();

  auto meshWithPointData =
      ls::SmartPointer<ls::Mesh<>>::New(); // Mesh with point data
  ls::ToMesh<double, D>(sphere2, meshWithPointData).apply();
  ls::VTKWriter<double>(meshWithPointData, "sphere2_LS_with_point_data.vtp")
      .apply();

  // Save mesh to file
  ls::VTKWriter<double>(mesh, "sparsefield.vtp").apply();

  // Get the calculated difference metrics
  double sumSquaredDifferences = compareSparseField.getSumSquaredDifferences();
  unsigned numPoints = compareSparseField.getNumPoints();
  double rmse = compareSparseField.getRMSE();

  std::cout << "Sphere 1 center: (" << origin1[0] << ", " << origin1[1] << ")"
            << std::endl;
  std::cout << "Sphere 2 center: (" << origin2[0] << ", " << origin2[1] << ")"
            << std::endl;
  std::cout << "Sum of squared differences: " << sumSquaredDifferences
            << std::endl;
  std::cout << "Number of points compared: " << numPoints << std::endl;
  std::cout << "RMSE: " << rmse << std::endl;

  // Compare to regular narrow band comparison for validation
  std::cout << "\nComparing with regular narrow band comparison:" << std::endl;

  // Use the original sphere objects directly instead of copying them
  ls::CompareNarrowBand<double, D> compareNarrowBand(sphere1, sphere2);
  compareNarrowBand.apply();

  std::cout << "Regular narrow band RMSE: " << compareNarrowBand.getRMSE()
            << std::endl;
  std::cout << "Regular narrow band points: "
            << compareNarrowBand.getNumPoints() << std::endl;

  // Test with range restrictions
  std::cout << "\nTesting with restricted ranges:" << std::endl;

  // Test with restricted X range
  compareSparseField.setOutputMesh(nullptr); // do not create mesh
  compareSparseField.clearXRange();
  compareSparseField.clearYRange();
  compareSparseField.setXRange(-5, 5);
  compareSparseField.apply();
  std::cout << "RMSE with X range [-5, 5]: " << compareSparseField.getRMSE()
            << std::endl;
  std::cout << "Number of points in X range: "
            << compareSparseField.getNumPoints() << std::endl;

  // // Test with restricted Y range
  // compareSparseField.clearXRange();
  // compareSparseField.setYRange(-5, 5);
  // compareSparseField.apply();
  // std::cout << "RMSE with Y range [-5, 5]: " << compareSparseField.getRMSE()
  //           << std::endl;
  // std::cout << "Number of points in Y range: "
  //           << compareSparseField.getNumPoints() << std::endl;

  // Test with both X and Y range restrictions
  compareSparseField.setXRange(-3, 3);
  compareSparseField.setYRange(-3, 3);
  compareSparseField.apply();
  std::cout << "RMSE with X range [-3, 3] and Y range [-3, 3]: "
            << compareSparseField.getRMSE() << std::endl;
  std::cout << "Number of points in both ranges: "
            << compareSparseField.getNumPoints() << std::endl;

  // Create a mesh output with squared differences
  compareSparseField.setOutputMesh(mesh);
  compareSparseField.apply();
  ls::VTKWriter<double>(mesh, "sparsefield_restricted.vtp").apply();

  // Clear range restrictions
  compareSparseField.clearXRange();
  compareSparseField.clearYRange();
  compareNarrowBand.clearXRange();
  compareNarrowBand.clearYRange();

  // Performance test: Compare time against regular narrow band comparison
  auto t1 = std::chrono::high_resolution_clock::now();

  compareSparseField.apply();

  auto t2 = std::chrono::high_resolution_clock::now();

  compareNarrowBand.apply();

  auto t3 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> sparse_ms = t2 - t1;
  std::chrono::duration<double, std::milli> narrowband_ms = t3 - t2;

  std::cout << "\nPerformance comparison:" << std::endl;
  std::cout << "Sparse Field execution time: " << sparse_ms.count() << " ms"
            << std::endl;
  std::cout << "Narrow Band execution time: " << narrowband_ms.count() << "ms "
            << std::endl;
  std::cout << "Performance ratio: "
            << narrowband_ms.count() / sparse_ms.count() << "x" << std::endl;

  return 0;
}
