#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsSliceExtractor.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsSliceExtractor that extracts a 2D slice from a 3D level set domain.
  This test creates a 3D sphere and extracts the z=0 plane as a 2D circle.
    \example SliceExtractor.cpp
*/

namespace ls = viennals;

int main() {
  // Create a 3D sphere
  double extent = 20;
  double gridDelta = 0.5;

  double bounds3D[6] = {-extent, extent, -extent, extent, -extent, extent};
  ls::Domain<double, 3>::BoundaryType boundaryCons3D[3];
  for (unsigned i = 0; i < 3; ++i)
    boundaryCons3D[i] =
        ls::Domain<double, 3>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere3D = ls::SmartPointer<ls::Domain<double, 3>>::New(
      bounds3D, boundaryCons3D, gridDelta);

  double origin[3] = {0., 0., 0.};
  double radius = 10.0;

  ls::MakeGeometry<double, 3>(
      sphere3D, ls::SmartPointer<ls::Sphere<double, 3>>::New(origin, radius))
      .apply();

  // Visualize 3D sphere
  auto mesh3D = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 3>(sphere3D, mesh3D).apply();
  ls::VTKWriter<double>(mesh3D, "sphere3D.vtp").apply();

  // Create a 2D domain to hold the slice
  double bounds2D[4] = {-extent, extent, -extent, extent};
  ls::Domain<double, 2>::BoundaryType boundaryCons2D[2];
  for (unsigned i = 0; i < 2; ++i)
    boundaryCons2D[i] =
        ls::Domain<double, 2>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto slice2D = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);

  // Extract the z=0 slice (dimension 2 = z-axis)
  ls::SliceExtractor<double> extractor;
  extractor.setSourceDomain(sphere3D);
  extractor.setSliceDomain(slice2D);
  extractor.setSliceDimension(2);  // z-axis
  extractor.setSlicePosition(0.0); // z=0 plane
  extractor.apply();

  auto mesh2D = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(slice2D, mesh2D).apply();
  ls::VTKWriter<double>(mesh2D, "slice2D.vtp").apply();

  // Test extracting different slices
  // Extract the x=5 slice (dimension 0 = x-axis)
  auto sliceX = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  ls::SliceExtractor<double>(sphere3D, sliceX, 0, 5.0).apply();
  auto meshX = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceX, meshX).apply();
  ls::VTKWriter<double>(meshX, "sliceX5.vtp").apply();

  // Extract the y=-5 slice (dimension 1 = y-axis)
  auto sliceY = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  ls::SliceExtractor<double>(sphere3D, sliceY, 1, -5.0).apply();
  auto meshY = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceY, meshY).apply();
  ls::VTKWriter<double>(meshY, "sliceY-5.vtp").apply();

  // Expand sliceX using lsExpand
  ls::Expand<double, 2>(sliceX, 10).apply();
  auto meshXExpanded = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceX, meshXExpanded).apply();
  ls::VTKWriter<double>(meshXExpanded, "sliceX5_expanded.vtp").apply();

  // Try slicing at a position that does not intersect the sphere
  auto sliceNoIntersection = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  ls::SliceExtractor<double>(sphere3D, sliceNoIntersection, 2, 15.0).apply();
  auto meshNoIntersection = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceNoIntersection, meshNoIntersection).apply();
  ls::VTKWriter<double>(meshNoIntersection, "sliceNoIntersection.vtp").apply();

  // Try with a large tolerance
  auto sliceLargeTolerance = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  auto extractorLT =
      ls::SliceExtractor<double>(sphere3D, sliceLargeTolerance, 2, 0.0);
  extractorLT.setTolerance(10.0);
  extractorLT.apply();
  auto meshLargeTolerance = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceLargeTolerance, meshLargeTolerance).apply();
  ls::VTKWriter<double>(meshLargeTolerance, "sliceLargeTolerance.vtp").apply();

  return 0;
}
