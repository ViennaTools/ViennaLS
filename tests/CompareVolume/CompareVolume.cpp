#include <cmath>
#include <iostream>
#include <string>

#include <lsCompareVolume.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsCompareArea that compares the area difference between two level
  sets. This test creates two different spheres and measures their area
  difference. \example CompareArea.cpp
*/

namespace ls = viennals;

template <int D> void runTest() {
  std::cout << "Running " << D << "D Test..." << std::endl;
  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D];
  for (int i = 0; i < 2 * D; ++i)
    bounds[i] = (i % 2 == 0) ? -extent : extent;

  typename ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  // Create first sphere (target)
  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  std::vector<double> origin1(D, 0.0);
  double radius1 = 5.0;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  // Create second sphere (sample) with different radius
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  std::vector<double> origin2(D, 0.0);
  double radius2 = 8.0;

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  std::string suffix = "_" + std::to_string(D) + "D.vtp";
  // Export both spheres as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere1" + suffix).apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere2, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere2" + suffix).apply();
  }

  // Compare the volumes/areas with mesh Output
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New(); // Create mesh for output

  if constexpr (D == 2) {
    ls::CompareArea<double> compareArea(sphere1, sphere2);
    compareArea.setOutputMesh(mesh);
    compareArea.apply();
    ls::VTKWriter<double>(mesh, "volumeDifference" + suffix + ".vtu").apply();

    // Area of circle = π * r²
    double area1 = M_PI * radius1 * radius1;
    double area2 = M_PI * radius2 * radius2;
    double theoreticalDiff = std::abs(area2 - area1);

    double calculatedDifference = compareArea.getAreaMismatch();
    unsigned long int cellCount = compareArea.getCellCount();

    std::cout << "Sphere 1 radius: " << radius1 << std::endl;
    std::cout << "Sphere 2 radius: " << radius2 << std::endl;
    std::cout << "Theoretical difference: " << theoreticalDiff << std::endl;
    std::cout << "Calculated difference: " << calculatedDifference << std::endl;
    std::cout << "Number of differing cells: " << cellCount << std::endl;
    std::cout << "Error: " << std::abs(calculatedDifference - theoreticalDiff)
              << std::endl;

    std::cout << "\nTesting custom increments and ranges:" << std::endl;
    compareArea.setDefaultIncrement(2);
    compareArea.apply();
    std::cout << "Difference with default increment of 2: "
              << compareArea.getCustomAreaMismatch() << std::endl;
    std::cout << "Cell count with default increment of 2: "
              << compareArea.getCustomCellCount() << std::endl;

    compareArea.setDefaultIncrement(1);
    compareArea.setXRangeAndIncrement(-5, 5, 3);
    compareArea.apply();
    std::cout << "Difference with x-range increment of 3: "
              << compareArea.getCustomAreaMismatch() << std::endl;
    std::cout << "Cell count with x-range increment of 3: "
              << compareArea.getCustomCellCount() << std::endl;

    compareArea.setDefaultIncrement(1);
    compareArea.setYRangeAndIncrement(-5, 5, 4);
    compareArea.apply();
    std::cout << "Difference with y-range increment of 4: "
              << compareArea.getCustomAreaMismatch() << std::endl;
    std::cout << "Cell count with y-range increment of 4: "
              << compareArea.getCustomCellCount() << std::endl;

  } else {
    ls::CompareVolume<double, D> compareVolume(sphere1, sphere2);
    compareVolume.setOutputMesh(mesh);
    compareVolume.apply();
    ls::VTKWriter<double>(mesh, "volumeDifference" + suffix + ".vtu").apply();

    // Volume of sphere = 4/3 * π * r³
    double vol1 = (4.0 / 3.0) * M_PI * std::pow(radius1, 3);
    double vol2 = (4.0 / 3.0) * M_PI * std::pow(radius2, 3);
    double theoreticalDiff = std::abs(vol2 - vol1);

    double calculatedDifference = compareVolume.getVolumeMismatch();
    unsigned long int cellCount = compareVolume.getCellCount();

    std::cout << "Sphere 1 radius: " << radius1 << std::endl;
    std::cout << "Sphere 2 radius: " << radius2 << std::endl;
    std::cout << "Theoretical difference: " << theoreticalDiff << std::endl;
    std::cout << "Calculated difference: " << calculatedDifference << std::endl;
    std::cout << "Number of differing cells: " << cellCount << std::endl;
    std::cout << "Error: " << std::abs(calculatedDifference - theoreticalDiff)
              << std::endl;

    std::cout << "\nTesting custom increments and ranges:" << std::endl;
    compareVolume.setDefaultIncrement(2);
    compareVolume.apply();
    std::cout << "Difference with default increment of 2: "
              << compareVolume.getCustomVolumeMismatch() << std::endl;
    std::cout << "Cell count with default increment of 2: "
              << compareVolume.getCustomCellCount() << std::endl;

    compareVolume.setDefaultIncrement(1);
    compareVolume.setXRangeAndIncrement(-5, 5, 3);
    compareVolume.apply();
    std::cout << "Difference with x-range increment of 3: "
              << compareVolume.getCustomVolumeMismatch() << std::endl;
    std::cout << "Cell count with x-range increment of 3: "
              << compareVolume.getCustomCellCount() << std::endl;

    compareVolume.setDefaultIncrement(1);
    compareVolume.setYRangeAndIncrement(-5, 5, 4);
    compareVolume.apply();
    std::cout << "Difference with y-range increment of 4: "
              << compareVolume.getCustomVolumeMismatch() << std::endl;
    std::cout << "Cell count with y-range increment of 4: "
              << compareVolume.getCustomCellCount() << std::endl;

    compareVolume.setDefaultIncrement(1);
    compareVolume.setZRangeAndIncrement(-5, 5, 5);
    compareVolume.apply();
    std::cout << "Difference with z-range increment of 5: "
              << compareVolume.getCustomVolumeMismatch() << std::endl;
    std::cout << "Cell count with z-range increment of 5: "
              << compareVolume.getCustomCellCount() << std::endl;
  }
}

int main() {
  omp_set_num_threads(4);
  runTest<2>();
  runTest<3>();
  return 0;
}
