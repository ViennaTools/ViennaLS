#include <iostream>

#include <lsCompareCriticalDimensions.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsCompareCriticalDimensions that compares critical dimensions
  (surface positions) between two level sets. This test creates two circles
  with shifted centers and compares their maximum and minimum positions.
  \example CompareCriticalDimensions.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  // Create first circle (reference)
  auto circle1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin1[D] = {0., 0.};
  double radius1 = 5.0;

  ls::MakeGeometry<double, D>(
      circle1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  // Create second circle (compare) with shifted center
  auto circle2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin2[D] = {1.5, 0.5}; // Shifted center
  double radius2 = 5.0;           // Same radius

  ls::MakeGeometry<double, D>(
      circle2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Export both circles as VTK files for visualization
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(circle1, mesh).apply();
    ls::VTKWriter<double>(mesh, "circle1_reference.vtp").apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(circle2, mesh).apply();
    ls::VTKWriter<double>(mesh, "circle2_compare.vtp").apply();
  }

  // Compare critical dimensions
  ls::CompareCriticalDimensions<double, D> compareCriticalDims(circle1,
                                                               circle2);

  // Add X ranges to find maximum and minimum Y positions (top and bottom)
  // Search in the central X range where both circles overlap
  compareCriticalDims.addXRange(-0, 0, true);  // Find maximum Y (top)
  compareCriticalDims.addXRange(-0, 0, false); // Find minimum Y (bottom)

  // Add Y ranges to find maximum and minimum X positions (right and left)
  // Search in the central Y range where both circles overlap
  compareCriticalDims.addYRange(-0, 0, true);  // Find maximum X (right)
  compareCriticalDims.addYRange(-0, 0, false); // Find minimum X (left)

  // Create mesh for output
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  compareCriticalDims.setOutputMesh(mesh);

  // Apply the comparison
  compareCriticalDims.apply();

  // Save mesh to file
  ls::VTKWriter<double>(mesh, "criticalDimensions.vtp").apply();

  // Debug: Print some surface mesh nodes to see actual positions
  std::cout << "\nDebug - Sample surface nodes from circle1:" << std::endl;
  auto debugMesh1 = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToSurfaceMesh<double, D>(circle1, debugMesh1).apply();
  std::cout << "Total nodes in circle1 surface: " << debugMesh1->nodes.size()
            << std::endl;
  for (size_t i = 0; i < std::min(size_t(10), debugMesh1->nodes.size()); ++i) {
    std::cout << "  Node " << i << ": (" << debugMesh1->nodes[i][0] << ", "
              << debugMesh1->nodes[i][1] << ")" << std::endl;
  }

  // Print results
  std::cout << "Circle 1 center: (" << origin1[0] << ", " << origin1[1] << ")"
            << std::endl;
  std::cout << "Circle 2 center: (" << origin2[0] << ", " << origin2[1] << ")"
            << std::endl;
  std::cout << "Radius: " << radius1 << std::endl;
  std::cout << "Center shift: (" << (origin2[0] - origin1[0]) << ", "
            << (origin2[1] - origin1[1]) << ")" << std::endl;
  std::cout << std::endl;

  // Get statistics
  size_t numCriticalDims = compareCriticalDims.getNumCriticalDimensions();
  std::cout << "Number of critical dimensions compared: " << numCriticalDims
            << std::endl;

  std::cout << "\nIndividual critical dimension results:" << std::endl;
  for (size_t i = 0; i < numCriticalDims; ++i) {
    double posRef, posCmp, diff;
    if (compareCriticalDims.getCriticalDimensionResult(i, posRef, posCmp,
                                                       diff)) {
      std::cout << "  Dimension " << i << ":" << std::endl;
      std::cout << "    Reference position: " << posRef << std::endl;
      std::cout << "    Compare position: " << posCmp << std::endl;
      std::cout << "    Difference: " << diff << std::endl;
    } else {
      std::cout << "  Dimension " << i << ": Invalid (not found)" << std::endl;
    }
  }

  std::cout << "\nAggregate statistics:" << std::endl;
  std::cout << "Mean difference: " << compareCriticalDims.getMeanDifference()
            << std::endl;
  std::cout << "Max difference: " << compareCriticalDims.getMaxDifference()
            << std::endl;
  std::cout << "RMSE: " << compareCriticalDims.getRMSE() << std::endl;

  // Theoretical validation
  // For circles with the same radius but shifted centers,
  // the critical dimensions should differ by approximately the shift amount
  // in the corresponding direction
  std::cout << "\nExpected differences based on geometry:" << std::endl;
  std::cout << "Top/Bottom (max/min Y): should be close to Y-shift = "
            << std::abs(origin2[1] - origin1[1]) << std::endl;
  std::cout << "Left/Right (max/min X): should be close to X-shift = "
            << std::abs(origin2[0] - origin1[0]) << std::endl;

  // Additional test: Test with wider ranges
  std::cout << "\n--- Testing with wider X range ---" << std::endl;
  compareCriticalDims.clearRanges();
  compareCriticalDims.addXRange(-10, 10, true);  // Find maximum Y
  compareCriticalDims.addXRange(-10, 10, false); // Find minimum Y
  compareCriticalDims.setOutputMesh(nullptr);    // Don't create mesh
  compareCriticalDims.apply();

  std::cout << "Number of critical dimensions: "
            << compareCriticalDims.getNumCriticalDimensions() << std::endl;
  for (size_t i = 0; i < compareCriticalDims.getNumCriticalDimensions(); ++i) {
    double posRef, posCmp, diff;
    if (compareCriticalDims.getCriticalDimensionResult(i, posRef, posCmp,
                                                       diff)) {
      std::cout << "  Dimension " << i << ": difference = " << diff
                << std::endl;
    }
  }

  // Additional test: Test with Y ranges only
  std::cout << "\n--- Testing with Y range only ---" << std::endl;
  compareCriticalDims.clearRanges();
  compareCriticalDims.addYRange(-10, 10, true);  // Find maximum X
  compareCriticalDims.addYRange(-10, 10, false); // Find minimum X
  compareCriticalDims.apply();

  std::cout << "Number of critical dimensions: "
            << compareCriticalDims.getNumCriticalDimensions() << std::endl;
  for (size_t i = 0; i < compareCriticalDims.getNumCriticalDimensions(); ++i) {
    double posRef, posCmp, diff;
    if (compareCriticalDims.getCriticalDimensionResult(i, posRef, posCmp,
                                                       diff)) {
      std::cout << "  Dimension " << i << ": difference = " << diff
                << std::endl;
    }
  }

  return 0;
}
