#include <cmath>
#include <iostream>

#include <lsCompareNarrowBand.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsCompareNarrowBand that compares two level sets by measuring
  the SDF value differences in their narrow bands. This test creates two
  spheres with different centers and measures their differences.
  \example CompareNarrowBand.cpp
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
  double radius2 = 5.0;         // Same radius

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Export both spheres as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere1_narrowband.vtp").apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere2, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere2_narrowband.vtp").apply();
  }

  // Compare the narrow bands
  ls::CompareNarrowBand<double, D> compareNarrowBand(sphere1, sphere2);

  // Create mesh for visualization of differences
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  compareNarrowBand.setOutputMesh(mesh, false);

  compareNarrowBand.apply();

  // Save mesh to file
  ls::VTKWriter<double>(mesh, "narrowband_absolute_differences.vtu").apply();

  // Get the calculated difference metrics
  double sumSquaredDifferences = compareNarrowBand.getSumSquaredDifferences();
  unsigned numPoints = compareNarrowBand.getNumPoints();
  double rmse = compareNarrowBand.getRMSE();

  std::cout << "Sphere 1 center: (" << origin1[0] << ", " << origin1[1] << ")"
            << std::endl;
  std::cout << "Sphere 2 center: (" << origin2[0] << ", " << origin2[1] << ")"
            << std::endl;
  std::cout << "Sum of squared differences: " << sumSquaredDifferences
            << std::endl;
  std::cout << "Number of points compared: " << numPoints << std::endl;
  std::cout << "RMSE: " << rmse << std::endl;

  // Test with range restrictions
  std::cout << "\nTesting with restricted ranges:" << std::endl;

  // Test with restricted X range
  compareNarrowBand.setOutputMesh(nullptr); // do not create mesh
  compareNarrowBand.clearXRange();
  compareNarrowBand.clearYRange();
  compareNarrowBand.setXRange(-5, 5);
  compareNarrowBand.apply();
  std::cout << "RMSE with X range [-5, 5]: " << compareNarrowBand.getRMSE()
            << std::endl;
  std::cout << "Number of points in X range: "
            << compareNarrowBand.getNumPoints() << std::endl;

  // Test with restricted Y range
  compareNarrowBand.clearXRange();
  compareNarrowBand.setYRange(-5, 5);
  compareNarrowBand.apply();
  std::cout << "RMSE with Y range [-5, 5]: " << compareNarrowBand.getRMSE()
            << std::endl;
  std::cout << "Number of points in Y range: "
            << compareNarrowBand.getNumPoints() << std::endl;

  // Test with both X and Y range restrictions
  compareNarrowBand.setXRange(-3, 3);
  compareNarrowBand.setYRange(-3, 3);
  compareNarrowBand.apply();
  std::cout << "RMSE with X range [-3, 3] and Y range [-3, 3]: "
            << compareNarrowBand.getRMSE() << std::endl;
  std::cout << "Number of points in both ranges: "
            << compareNarrowBand.getNumPoints() << std::endl;

  // Create a mesh output with squared differences
  compareNarrowBand.setOutputMesh(mesh);
  compareNarrowBand.setOutputMeshSquaredDifferences(true);
  compareNarrowBand.apply();
  ls::VTKWriter<double>(mesh,
                        "narrowband_resctricted-range_squared_differences.vtu")
      .apply();

  return 0;
}
