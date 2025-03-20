#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  emulation, by creating a trench geometry. A uniform
  layer of a different material is then grown on top. It is
  the same example as Deposition but emulates the deposition
  rather than simulating a slow growth.
  \example GeometricAdvection.cpp
*/

namespace ls = viennals;

using NumericType = float;

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  NumericType extent = 30;
  NumericType gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  {
    NumericType origin[3] = {0., 0., 0.};
    NumericType planeNormal[3] = {0., 0., 1.};
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    // make -x and +x greater than domain for numerical stability
    NumericType ylimit = extent / 4.;
    NumericType minCorner[D] = {-extent - 1, -ylimit, -15.};
    NumericType maxCorner[D] = {extent + 1, ylimit, 1.};
    auto box =
        ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<NumericType, D>(trench, box).apply();
    // Create trench geometry
    ls::BooleanOperation<NumericType, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "trench-0.vtp").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  auto newLayer = ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);

  std::cout << "Advecting" << std::endl;
  // Grow the layer uniformly by 4 as in deposition example
  auto dist =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          4.0, gridDelta);
  ls::GeometricAdvect<NumericType, D>(newLayer, dist).apply();

  {
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "trench-final.vtp").apply();
  }

  return 0;
}
