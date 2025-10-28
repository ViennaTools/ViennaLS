#include <chrono>
#include <iostream>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsCalculateCurvatures.hpp>

#include <omp.h>

/**
  Minimal example of how to calculate the Curvatures of the level set function.
  Works in 2D and 3D.
  Outputs the selected Curvatures into a .vtk file.
  \example calculateCurvatures.cpp
*/

namespace ls = viennals;

constexpr int D = 3;
typedef double NumericType;

int main() {

  omp_set_num_threads(1);

  NumericType gridDelta = 0.5;

  auto sphere = ls::SmartPointer<ls::Domain<NumericType, D>>::New(gridDelta);
  NumericType origin[3] = {5., 0., 0.};
  NumericType radius = 10.0;

  ls::MakeGeometry<NumericType, D>(
      sphere, ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
      .apply();

  ls::Expand<NumericType, D>(sphere, 5).apply();

  ls::CalculateCurvatures<NumericType, D> calcCurve(sphere);

  if (D == 2) {
    calcCurve.setCurvatureType(ls::CurvatureEnum::MEAN_CURVATURE);
  } else {
    calcCurve.setCurvatureType(ls::CurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE);
  }

  calcCurve.apply();

  auto meanCurvatures = sphere->getPointData().getScalarData("MeanCurvatures");

  VC_TEST_ASSERT(meanCurvatures != nullptr)

  double analyticCurvature = 1. / radius;
  viennahrle::SizeType numberOfActivePoints = 0;
  double sum = 0.;
  for (viennahrle::ConstSparseIterator<
           typename ls::Domain<NumericType, D>::DomainType>
           it(sphere->getDomain());
       !it.isFinished(); ++it) {
    if (NumericType value = it.getValue();
        !it.isDefined() || std::abs(value) > 0.5)
      continue;

    sum += meanCurvatures->at(it.getPointId());
    ++numberOfActivePoints;
  }

  VC_TEST_ASSERT(std::abs(sum / numberOfActivePoints - analyticCurvature) <
                 1e-3)

  // std::cout << "Writing Output..." << std::endl;
  // typename lsPointData<NumericType>::ScalarDataType analytic;
  // analytic.resize(meanCurvatures->size(), analyticCurvature);
  // sphere->getPointData().insertNextScalarData(analytic,
  // "analyticCurvatures");

  // auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  // lsToSurfaceMesh(sphere, mesh).apply();
  // ls::VTKWriter(mesh, "curvatures.vtp").apply();

  // ls::ToMesh<NumericType, D>(sphere, mesh, true, true).apply();
  // auto writer = ls::VTKWriter<NumericType>();
  // writer.setMesh(mesh);
  // writer.setFileName("curvatures.vtk");
  // writer.apply();

  // std::cout << "Finished" << std::endl;

  return 0;
}
