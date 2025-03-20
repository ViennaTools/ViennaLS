#include <iostream>

#include <lsCompareArea.hpp>
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

  // Create second sphere (sample) with different radius
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin2[D] = {0., 0.};
  double radius2 = 8.0;

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Export both spheres as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere1.vtp").apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere2, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere2.vtp").apply();
  }

  // Compare the areas with mesh Output
  ls::CompareArea<double, D> compareArea(sphere1, sphere2);
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New(); // Create mesh for output
  compareArea.setOutputMesh(mesh);
  compareArea.apply();
  // save mesh to file
  ls::VTKWriter<double>(mesh, "areaDifference.vtu").apply();

  // Calculate theoretical area difference (for circles in 2D)
  // Area of circle = π * r²
  double theoreticalArea1 = M_PI * radius1 * radius1;
  double theoreticalArea2 = M_PI * radius2 * radius2;
  double theoreticalDifference = std::abs(theoreticalArea2 - theoreticalArea1);

  // Get the calculated area difference
  double calculatedDifference = compareArea.getAreaMismatch();
  unsigned long int cellCount = compareArea.getCellCount();

  std::cout << "Sphere 1 radius: " << radius1 << std::endl;
  std::cout << "Sphere 2 radius: " << radius2 << std::endl;
  std::cout << "Theoretical area difference: " << theoreticalDifference
            << std::endl;
  std::cout << "Calculated area difference: " << calculatedDifference
            << std::endl;
  std::cout << "Number of differing cells: " << cellCount << std::endl;
  std::cout << "Error: "
            << std::abs(calculatedDifference - theoreticalDifference)
            << std::endl;

  // Test custom increment and range functionality
  std::cout << "\nTesting custom increments and ranges:" << std::endl;

  // Set custom increment for whole domain
  compareArea.setDefaultIncrement(2);
  compareArea.apply();
  std::cout << "Area difference with default increment of 2: "
            << compareArea.getCustomAreaMismatch() << std::endl;
  std::cout << "Cell count with default increment of 2: "
            << compareArea.getCustomCellCount() << std::endl;

  // Set range-specific increment for x-range
  compareArea.setDefaultIncrement(1);
  compareArea.setXRangeAndIncrement(-5, 5, 3);
  compareArea.apply();
  std::cout << "Area difference with x-range increment of 3: "
            << compareArea.getCustomAreaMismatch() << std::endl;
  std::cout << "Cell count with x-range increment of 3: "
            << compareArea.getCustomCellCount() << std::endl;

  // Set range-specific increment for y-range
  compareArea.setDefaultIncrement(1);
  compareArea.setYRangeAndIncrement(-5, 5, 4);
  compareArea.apply();
  std::cout << "Area difference with y-range increment of 4: "
            << compareArea.getCustomAreaMismatch() << std::endl;
  std::cout << "Cell count with y-range increment of 4: "
            << compareArea.getCustomCellCount() << std::endl;

  return 0;
}
