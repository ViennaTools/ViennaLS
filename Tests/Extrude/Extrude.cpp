#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExtrude.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Simple example showing how to use lsExtrude in order
  to transform a 2D trench into a 3D level set.
  \example Extrude.cpp
*/

int main() {

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  // 2D domain boundaries
  double bounds[2 * 2] = {-extent, extent, -extent, extent};
  lsDomain<double, 2>::BoundaryType boundaryCons[2];
  boundaryCons[0] = lsDomain<double, 2>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<double, 2>::BoundaryType::INFINITE_BOUNDARY;

  auto trench =
      lsSmartPointer<lsDomain<double, 2>>::New(bounds, boundaryCons, gridDelta);

  {
    double origin[2] = {0., 0.};
    double normal[2] = {0., 1.};

    lsMakeGeometry<double, 2>(
        trench, lsSmartPointer<lsPlane<double, 2>>::New(origin, normal))
        .apply();

    auto cutOut = lsSmartPointer<lsDomain<double, 2>>::New(bounds, boundaryCons,
                                                           gridDelta);

    double minPoint[2] = {-5., -5.};
    double maxPoint[2] = {5., gridDelta};

    lsMakeGeometry<double, 2>(
        cutOut, lsSmartPointer<lsBox<double, 2>>::New(minPoint, maxPoint))
        .apply();

    lsBooleanOperation<double, 2>(trench, cutOut,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<double, 2>(trench, mesh).apply();
    lsVTKWriter<double>(mesh, "trench_initial.vtp").apply();
  }

  std::array<double, 2> extrudeExtent = {-5., 5.};
  auto trench_3D = lsSmartPointer<lsDomain<double, 3>>::New();
  lsExtrude<double>(trench, trench_3D, extrudeExtent, 1).apply();

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<double, 3>(trench_3D, mesh).apply();
    lsVTKWriter<double>(mesh, "trench_extrude.vtp").apply();
  }

  return 0;
}
