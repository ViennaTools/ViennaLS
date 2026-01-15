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

template <int D> void runTest() {
  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D];
  for (unsigned i = 0; i < D; ++i) {
    bounds[2 * i] = -extent;
    bounds[2 * i + 1] = extent;
  }

  typename ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  // Create first sphere (target)
  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin1[D];
  for (int i = 0; i < D; ++i)
    origin1[i] = 0.;
  double radius1 = 5.0;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  // Create second sphere (sample) with different center
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin2[D];
  for (int i = 0; i < D; ++i)
    origin2[i] = 0.;
  origin2[0] = 2.;
  if (D > 1)
    origin2[1] = 1.;
  // double origin2[D] = {0., 0.}; // Same center
  double radius2 = 5.0; // Same radius

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Make sure target is expanded
  // Note: We'll let CompareSparseField do the expansion automatically
  // ls::Expand<double, D>(sphere1, 50).apply();

  // Reduce the sample level set to a sparse field
  ls::Reduce<double, D>(sphere2, 1).apply();

  std::string dimString = std::to_string(D) + "D";

  // Export both spheres as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere1_expanded_" + dimString + ".vtp")
        .apply();
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere1, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere1_surface_" + dimString + ".vtp")
        .apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere2, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere2_sparse_iterated_" + dimString + ".vtp")
        .apply();
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere2, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere2_surface_" + dimString + ".vtp")
        .apply();
  }

  // Compare using sparse field comparison
  ls::CompareSparseField<double, D> compareSparseField(sphere1, sphere2);

  // Test the new setExpandedLevelSetWidth feature
  // Set a custom expansion width (default is 50)
  compareSparseField.setExpandedLevelSetWidth(75);
  std::cout << "Using custom expansion width of 75 for the expanded level set"
            << std::endl;

  // Create mesh for visualization of differences
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  compareSparseField.setFillIteratedWithDistances(true);
  compareSparseField.setOutputMesh(mesh);
  compareSparseField.apply();

  auto meshWithPointData =
      ls::SmartPointer<ls::Mesh<>>::New(); // Mesh with point data
  ls::ToMesh<double, D>(sphere2, meshWithPointData).apply();
  ls::VTKWriter<double>(meshWithPointData,
                        "sphere2_LS_with_point_data_" + dimString + ".vtp")
      .apply();

  // Save mesh to file
  ls::VTKWriter<double>(mesh, "sparsefield_" + dimString + ".vtp").apply();

  // Get the calculated difference metrics
  double sumSquaredDifferences = compareSparseField.getSumSquaredDifferences();
  double rmse = compareSparseField.getRMSE();

  // Check number of points
  unsigned numPoints = compareSparseField.getNumPoints();
  unsigned numSkippedPoints =
      compareSparseField.getNumSkippedPoints(); // Number of skipped points

  std::cout << "\nComparison Results (" << dimString << "):" << std::endl;
  std::cout << "Sphere 1 center: (";
  for (int i = 0; i < D; ++i)
    std::cout << origin1[i] << ((i == D - 1) ? "" : ", ");
  std::cout << ")" << std::endl;
  std::cout << "Sphere 2 center: (";
  for (int i = 0; i < D; ++i)
    std::cout << origin2[i] << ((i == D - 1) ? "" : ", ");
  std::cout << ")" << std::endl;
  std::cout << "Sphere 1 level set width after expansion: "
            << sphere1->getLevelSetWidth() << std::endl;
  std::cout << "Sum of squared differences: " << sumSquaredDifferences
            << std::endl;
  std::cout << "Number of points compared: " << numPoints << std::endl;
  std::cout << "RMSE: " << rmse << std::endl;
  std::cout << "Number of skipped points: " << numSkippedPoints << std::endl;

  // Compare to regular narrow band comparison for validation
  std::cout << "\nComparing with regular narrow band comparison:" << std::endl;

  // // Use the original sphere objects directly instead of copying them
  // ls::CompareNarrowBand<double, D> compareNarrowBand(sphere1, sphere2);
  // compareNarrowBand.apply();

  // std::cout << "Regular narrow band RMSE: " << compareNarrowBand.getRMSE()
  //           << std::endl;
  // std::cout << "Regular narrow band points: "
  //           << compareNarrowBand.getNumPoints() << std::endl;

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

  // Test with restricted Y range
  compareSparseField.clearXRange();
  compareSparseField.setYRange(-5, 5);
  compareSparseField.apply();
  std::cout << "RMSE with Y range [-5, 5]: " << compareSparseField.getRMSE()
            << std::endl;
  std::cout << "Number of points in Y range: "
            << compareSparseField.getNumPoints() << std::endl;

  if constexpr (D == 3) {
    // Test with restricted Z range
    compareSparseField.clearYRange();
    compareSparseField.setZRange(-5, 5);
    compareSparseField.apply();
    std::cout << "RMSE with Z range [-5, 5]: " << compareSparseField.getRMSE()
              << std::endl;
    std::cout << "Number of points in Z range: "
              << compareSparseField.getNumPoints() << std::endl;
    compareSparseField.clearZRange();
  }

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
  ls::VTKWriter<double>(mesh, "sparsefield_restricted_" + dimString + ".vtp")
      .apply();

  // Test with different expansion widths
  std::cout << "\nTesting with different expansion widths:" << std::endl;

  // Reset ranges
  compareSparseField.clearXRange();
  compareSparseField.clearYRange();
  compareSparseField.setOutputMesh(nullptr);

  // Test with smaller expansion width (should still work but may have less
  // coverage)
  auto sphere1_narrow = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<double, D>(
      sphere1_narrow,
      ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  ls::CompareSparseField<double, D> compareSmallWidth(sphere1_narrow, sphere2);
  compareSmallWidth.setExpandedLevelSetWidth(30);
  compareSmallWidth.apply();
  std::cout << "RMSE with expansion width 30: " << compareSmallWidth.getRMSE()
            << std::endl;
  std::cout << "Level set width after expansion: "
            << sphere1_narrow->getLevelSetWidth() << std::endl;

  // Test with larger expansion width
  auto sphere1_wide = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<double, D>(
      sphere1_wide,
      ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  ls::CompareSparseField<double, D> compareLargeWidth(sphere1_wide, sphere2);
  compareLargeWidth.setExpandedLevelSetWidth(100);
  compareLargeWidth.apply();
  std::cout << "RMSE with expansion width 100: " << compareLargeWidth.getRMSE()
            << std::endl;
  std::cout << "Level set width after expansion: "
            << sphere1_wide->getLevelSetWidth() << std::endl;

  // // Clear range restrictions
  // compareSparseField.clearXRange();
  // compareSparseField.clearYRange();
  // compareNarrowBand.clearXRange();
  // compareNarrowBand.clearYRange();

  // // Performance test: Compare time against regular narrow band comparison
  // auto t1 = std::chrono::high_resolution_clock::now();

  // compareSparseField.apply();

  // auto t2 = std::chrono::high_resolution_clock::now();

  // compareNarrowBand.apply();

  // auto t3 = std::chrono::high_resolution_clock::now();

  // std::chrono::duration<double, std::milli> sparse_ms = t2 - t1;
  // std::chrono::duration<double, std::milli> narrowband_ms = t3 - t2;

  // std::cout << "\nPerformance comparison:" << std::endl;
  // std::cout << "Sparse Field execution time: " << sparse_ms.count() << " ms"
  //           << std::endl;
  // std::cout << "Narrow Band execution time: " << narrowband_ms.count() << "ms
  // "
  //           << std::endl;
  // std::cout << "Performance ratio: "
  //           << narrowband_ms.count() / sparse_ms.count() << "x" << std::endl;
}

int main() {
  omp_set_num_threads(8);
  runTest<2>();
  runTest<3>();
  return 0;
}
