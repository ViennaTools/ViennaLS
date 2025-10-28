#include <chrono>
#include <cmath>
#include <iostream>

#include <lsCompareArea.hpp>
#include <lsCompareChamfer.hpp>
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
  Test for lsCompareChamfer that compares two level sets by computing the
  Chamfer distance between their surfaces. This test creates two circles with
  different centers and measures their geometric similarity using Chamfer
  distance, comparing it with other metrics.
  \example CompareChamfer.cpp
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

  // Create first circle (target)
  auto circle1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin1[D] = {0., 0.};
  double radius1 = 5.0;

  ls::MakeGeometry<double, D>(
      circle1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  // Create second circle (sample) with different center
  auto circle2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin2[D] = {2., 1.}; // Shifted center
  // double origin2[D] = {0., 0.}; // Same center (for testing)
  double radius2 = 5.0; // Same radius

  ls::MakeGeometry<double, D>(
      circle2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Export both circles as VTK files for visualization
  std::cout << "Exporting surface meshes..." << std::endl;
  {
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(circle1, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "circle1_surface.vtp").apply();
    std::cout << "  Circle 1 surface points: " << meshSurface->nodes.size()
              << std::endl;
  }

  {
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(circle2, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "circle2_surface.vtp").apply();
    std::cout << "  Circle 2 surface points: " << meshSurface->nodes.size()
              << std::endl;
  }

  // Test 1: Basic Chamfer distance calculation
  std::cout << "\n=== Test 1: Basic Chamfer Distance ===" << std::endl;
  std::cout << "Circle 1 center: (" << origin1[0] << ", " << origin1[1] << ")"
            << std::endl;
  std::cout << "Circle 2 center: (" << origin2[0] << ", " << origin2[1] << ")"
            << std::endl;
  std::cout << "Expected geometric shift: "
            << std::sqrt((origin2[0] - origin1[0]) * (origin2[0] - origin1[0]) +
                         (origin2[1] - origin1[1]) * (origin2[1] - origin1[1]))
            << std::endl;

  ls::CompareChamfer<double, D> compareChamfer(circle1, circle2);

  // Create output meshes with distance information
  auto targetMesh = ls::SmartPointer<ls::Mesh<>>::New();
  auto sampleMesh = ls::SmartPointer<ls::Mesh<>>::New();
  compareChamfer.setOutputMeshTarget(targetMesh);
  compareChamfer.setOutputMeshSample(sampleMesh);

  auto t1 = std::chrono::high_resolution_clock::now();
  compareChamfer.apply();
  auto t2 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> chamfer_ms = t2 - t1;

  std::cout << "\nChamfer Distance Results:" << std::endl;
  std::cout << "  Forward distance (target → sample): "
            << compareChamfer.getForwardDistance() << std::endl;
  std::cout << "  Backward distance (sample → target): "
            << compareChamfer.getBackwardDistance() << std::endl;
  std::cout << "  Chamfer distance (average): "
            << compareChamfer.getChamferDistance() << std::endl;
  std::cout << "  RMS Chamfer distance: "
            << compareChamfer.getRMSChamferDistance() << std::endl;
  std::cout << "  Maximum distance: " << compareChamfer.getMaxDistance()
            << std::endl;
  std::cout << "  Target surface points: "
            << compareChamfer.getNumTargetPoints() << std::endl;
  std::cout << "  Sample surface points: "
            << compareChamfer.getNumSamplePoints() << std::endl;
  std::cout << "  Execution time: " << chamfer_ms.count() << " ms" << std::endl;

  // Save meshes with distance data
  ls::VTKWriter<double>(targetMesh, "chamfer_target_distances.vtp").apply();
  ls::VTKWriter<double>(sampleMesh, "chamfer_sample_distances.vtp").apply();

  // Test 2: Compare with other metrics
  std::cout << "\n=== Test 2: Comparison with Other Metrics ===" << std::endl;

  // Sparse Field comparison
  auto circle1_expanded = ls::SmartPointer<ls::Domain<double, D>>::New(circle1);
  ls::Expand<double, D>(circle1_expanded, 50).apply();
  auto circle2_reduced = ls::SmartPointer<ls::Domain<double, D>>::New(circle2);
  ls::Reduce<double, D>(circle2_reduced, 1).apply();

  ls::CompareSparseField<double, D> compareSparseField(circle1_expanded,
                                                       circle2_reduced);
  auto t3 = std::chrono::high_resolution_clock::now();
  compareSparseField.apply();
  auto t4 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> sparse_ms = t4 - t3;

  std::cout << "Sparse Field Results:" << std::endl;
  std::cout << "  RMSE: " << compareSparseField.getRMSE() << std::endl;
  std::cout << "  Points compared: " << compareSparseField.getNumPoints()
            << std::endl;
  std::cout << "  Execution time: " << sparse_ms.count() << " ms" << std::endl;

  // Area comparison
  ls::CompareArea<double, D> compareArea(circle1, circle2);
  auto t5 = std::chrono::high_resolution_clock::now();
  compareArea.apply();
  auto t6 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> area_ms = t6 - t5;

  std::cout << "\nArea Comparison Results:" << std::endl;
  std::cout << "  Area mismatch: " << compareArea.getAreaMismatch()
            << std::endl;
  std::cout << "  Different cells: " << compareArea.getCellCount() << std::endl;
  std::cout << "  Execution time: " << area_ms.count() << " ms" << std::endl;

  // Test 3: Different geometric configurations
  std::cout << "\n=== Test 3: Different Geometric Configurations ==="
            << std::endl;

  // Test 3a: Identical circles (should give near-zero Chamfer distance)
  auto circle3 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<double, D>(
      circle3, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  ls::CompareChamfer<double, D> compareIdentical(circle1, circle3);
  compareIdentical.apply();

  std::cout << "Identical circles:" << std::endl;
  std::cout << "  Chamfer distance: " << compareIdentical.getChamferDistance()
            << " (expected ~0)" << std::endl;

  // Test 3b: Different radii
  auto circle4 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  double radius4 = 7.0; // Larger radius
  ls::MakeGeometry<double, D>(
      circle4, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius4))
      .apply();

  ls::CompareChamfer<double, D> compareDifferentSize(circle1, circle4);
  compareDifferentSize.apply();

  std::cout << "\nDifferent radii (r1=" << radius1 << ", r2=" << radius4
            << "):" << std::endl;
  std::cout << "  Chamfer distance: "
            << compareDifferentSize.getChamferDistance() << std::endl;
  std::cout << "  Expected difference: " << std::abs(radius4 - radius1)
            << std::endl;
  std::cout << "  Forward distance: "
            << compareDifferentSize.getForwardDistance() << std::endl;
  std::cout << "  Backward distance: "
            << compareDifferentSize.getBackwardDistance() << std::endl;

  // Test 3c: Large shift
  auto circle5 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  double origin5[D] = {5., 0.}; // Larger shift
  ls::MakeGeometry<double, D>(
      circle5, ls::SmartPointer<ls::Sphere<double, D>>::New(origin5, radius1))
      .apply();

  ls::CompareChamfer<double, D> compareLargeShift(circle1, circle5);
  compareLargeShift.apply();

  std::cout << "\nLarge shift (5 units in x-direction):" << std::endl;
  std::cout << "  Chamfer distance: " << compareLargeShift.getChamferDistance()
            << std::endl;
  std::cout << "  Expected shift: " << origin5[0] - origin1[0] << std::endl;

  // Test 4: Performance summary
  std::cout << "\n=== Performance Summary ===" << std::endl;
  std::cout << "Chamfer distance: " << chamfer_ms.count() << " ms" << std::endl;
  std::cout << "Sparse field:     " << sparse_ms.count() << " ms" << std::endl;
  std::cout << "Area comparison:  " << area_ms.count() << " ms" << std::endl;

  std::cout << "\n=== All Tests Completed Successfully ===" << std::endl;

  return 0;
}
