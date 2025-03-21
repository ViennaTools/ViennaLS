#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsSlice.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Test for lsSlice that extracts a 2D slice from a 3D level set domain.
  This test creates a 3D sphere and extracts the z=0 plane as a 2D circle.
    \example Slice.cpp
*/

namespace ls = viennals;

int main() {
  // Create a 3D sphere
  constexpr double extent = 20;
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
  ls::Slice<double> extractor;
  extractor.setSourceLevelSet(sphere3D);
  extractor.setSliceLevelSet(slice2D);
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
  ls::Slice<double>(sphere3D, sliceX, 0, 5.0).apply();
  auto meshX = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceX, meshX).apply();
  ls::VTKWriter<double>(meshX, "sliceX5.vtp").apply();

  // Extract the y=-5 slice (dimension 1 = y-axis)
  auto sliceY = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  ls::Slice<double>(sphere3D, sliceY, 1, -5.0).apply();
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
  ls::Slice<double>(sphere3D, sliceNoIntersection, 2, 15.0).apply();
  auto meshNoIntersection = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceNoIntersection, meshNoIntersection).apply();
  ls::VTKWriter<double>(meshNoIntersection, "sliceNoIntersection.vtp").apply();

  // Try slicing at a position not divisible by grid delta
  auto sliceNotDivisible = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2D, boundaryCons2D, gridDelta);
  ls::Slice<double>(sphere3D, sliceNotDivisible, 2, 2.75).apply();
  auto meshNotDivisible = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(sliceNotDivisible, meshNotDivisible).apply();
  ls::VTKWriter<double>(meshNotDivisible, "sliceNotDivisible.vtp").apply();

  // Try not passing a slice domain
  ls::Slice<double> extractorNoSliceDomain;
  extractorNoSliceDomain.setSourceLevelSet(sphere3D);
  extractorNoSliceDomain.setSliceDimension(2);  // z-axis
  extractorNoSliceDomain.setSlicePosition(0.0); // z=0 plane
  extractorNoSliceDomain.apply();
  auto meshNoSliceDomain = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(extractorNoSliceDomain.getSliceLevelSet(),
                        meshNoSliceDomain)
      .apply();
  ls::VTKWriter<double>(meshNoSliceDomain, "sliceNoSliceDomain.vtp").apply();

  // Try with bounds which end at 0
  double bounds2DZero[4] = {-extent, 0, -extent, 0};
  auto slice2DZero = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds2DZero, boundaryCons2D, gridDelta);
  ls::Slice<double>(sphere3D, slice2DZero, 2, 0.0).apply();
  auto mesh2DZero = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, 2>(slice2DZero, mesh2DZero).apply();
  ls::VTKWriter<double>(mesh2DZero, "slice2DZero.vtp").apply();

  return 0;
}
