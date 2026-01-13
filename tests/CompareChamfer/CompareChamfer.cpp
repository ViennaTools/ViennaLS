#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

#include <lsCompareChamfer.hpp>
#include <lsCompareSparseField.hpp>
#include <lsCompareVolume.hpp>
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

  // Create second sphere (sample) with different center
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  std::vector<double> origin2(D, 0.0);
  origin2[0] = 2.;
  origin2[1] = 1.;
  double radius2 = 5.0; // Same radius

  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin2, radius2))
      .apply();

  // Export both spheres as VTK files for visualization
  std::string suffix = "_" + std::to_string(D) + "D.vtp";
  std::cout << "Exporting surface meshes to *" << suffix << "..." << std::endl;
  {
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere1, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere1_surface" + suffix).apply();
    std::cout << "  Sphere 1 surface points: " << meshSurface->nodes.size()
              << std::endl;
  }

  {
    auto meshSurface = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere2, meshSurface).apply();
    ls::VTKWriter<double>(meshSurface, "sphere2_surface" + suffix).apply();
    std::cout << "  Sphere 2 surface points: " << meshSurface->nodes.size()
              << std::endl;
  }

  // Test 1: Basic Chamfer distance calculation
  std::cout << "\n=== Test 1: Basic Chamfer Distance ===" << std::endl;
  std::cout << "Sphere 1 center: (" << origin1[0] << ", " << origin1[1] << ")"
            << std::endl;
  std::cout << "Sphere 2 center: (" << origin2[0] << ", " << origin2[1] << ")"
            << std::endl;
  double distSq = 0.0;
  for (int i = 0; i < D; ++i)
    distSq += (origin2[i] - origin1[i]) * (origin2[i] - origin1[i]);
  std::cout << "Expected geometric shift: " << std::sqrt(distSq) << std::endl;

  ls::CompareChamfer<double, D> compareChamfer(sphere1, sphere2);

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
  ls::VTKWriter<double>(targetMesh, "chamfer_target_distances" + suffix)
      .apply();
  ls::VTKWriter<double>(sampleMesh, "chamfer_sample_distances" + suffix)
      .apply();

  // Test 2: Compare with other metrics
  std::cout << "\n=== Test 2: Comparison with Other Metrics ===" << std::endl;

  // Sparse Field comparison
  auto sphere1_expanded = ls::SmartPointer<ls::Domain<double, D>>::New(sphere1);
  ls::Expand<double, D>(sphere1_expanded, 50).apply();
  auto sphere2_reduced = ls::SmartPointer<ls::Domain<double, D>>::New(sphere2);
  ls::Reduce<double, D>(sphere2_reduced, 1).apply();

  ls::CompareSparseField<double, D> compareSparseField(sphere1_expanded,
                                                       sphere2_reduced);
  auto t3 = std::chrono::high_resolution_clock::now();
  compareSparseField.apply();
  auto t4 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> sparse_ms = t4 - t3;

  std::cout << "Sparse Field Results:" << std::endl;
  std::cout << "  RMSE: " << compareSparseField.getRMSE() << std::endl;
  std::cout << "  Points compared: " << compareSparseField.getNumPoints()
            << std::endl;
  std::cout << "  Execution time: " << sparse_ms.count() << " ms" << std::endl;

  // Area/Volume comparison
  ls::CompareVolume<double, D> compareVolume(sphere1, sphere2);
  auto t5 = std::chrono::high_resolution_clock::now();
  compareVolume.apply();
  auto t6 = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double, std::milli> area_ms = t6 - t5;

  std::cout << "\nArea/Volume Comparison Results:" << std::endl;
  std::cout << "  Area/Volume mismatch: " << compareVolume.getVolumeMismatch()
            << std::endl;
  std::cout << "  Different cells: " << compareVolume.getCellCount()
            << std::endl;
  std::cout << "  Execution time: " << area_ms.count() << " ms" << std::endl;

  // Test 3: Different geometric configurations
  std::cout << "\n=== Test 3: Different Geometric Configurations ==="
            << std::endl;

  // Test 3a: Identical spheres (should give near-zero Chamfer distance)
  auto sphere3 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<double, D>(
      sphere3, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius1))
      .apply();

  ls::CompareChamfer<double, D> compareIdentical(sphere1, sphere3);
  compareIdentical.apply();

  std::cout << "Identical spheres:" << std::endl;
  std::cout << "  Chamfer distance: " << compareIdentical.getChamferDistance()
            << " (expected ~0)" << std::endl;

  // Test 3b: Different radii
  auto sphere4 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  double radius4 = 7.0; // Larger radius
  ls::MakeGeometry<double, D>(
      sphere4, ls::SmartPointer<ls::Sphere<double, D>>::New(origin1, radius4))
      .apply();

  ls::CompareChamfer<double, D> compareDifferentSize(sphere1, sphere4);
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
  auto sphere5 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  std::vector<double> origin5(D, 0.0);
  origin5[0] = 5.; // Larger shift
  ls::MakeGeometry<double, D>(
      sphere5, ls::SmartPointer<ls::Sphere<double, D>>::New(origin5, radius1))
      .apply();

  ls::CompareChamfer<double, D> compareLargeShift(sphere1, sphere5);
  compareLargeShift.apply();

  std::cout << "\nLarge shift (5 units in x-direction):" << std::endl;
  std::cout << "  Chamfer distance: " << compareLargeShift.getChamferDistance()
            << std::endl;
  std::cout << "  Expected shift: " << origin5[0] - origin1[0] << std::endl;

  // Test 4: Performance summary
  std::cout << "\n=== Performance Summary ===" << std::endl;
  std::cout << "Chamfer distance: " << chamfer_ms.count() << " ms" << std::endl;
}

int main() {
  omp_set_num_threads(8);

  runTest<2>();
  runTest<3>();

  return 0;
}
