#include <chrono>
#include <iostream>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsCalculateCurvatures.hpp>

#include <omp.h>

/**
  Minimal example of how to calculate the Curvatures of the level set function.
  Works in 2D and 3D.
  Outputs the selected Curvatures into a .vtk file.
  \example calculateCurvatures.cpp
*/

constexpr int D = 3;
typedef double NumericType;

int main() {

  omp_set_num_threads(1);

  NumericType gridDelta = 0.5;

  std::cout << "Creating Sphere..." << std::endl;
  auto sphere = lsSmartPointer<lsDomain<NumericType, D>>::New(gridDelta);
  NumericType origin[3] = {5., 0., 0.};
  NumericType radius = 10.0;

  lsMakeGeometry<NumericType, D>(
      sphere, lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();

  std::cout << "Expanding..." << std::endl;
  lsExpand<NumericType, D>(sphere, 5).apply();

  std::cout << "Calculating Curvatures..." << std::endl;

  lsCalculateCurvatures<NumericType, D> calcCurve(sphere);

  if (D == 2) {
    calcCurve.setCurvatureType(lsCurvatureType::MEAN_CURVATURE);
  } else {
    calcCurve.setCurvatureType(lsCurvatureType::MEAN_AND_GAUSSIAN_CURVATURE);
  }

  calcCurve.apply();

  std::cout << "Writing Output..." << std::endl;

  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  lsToMesh<NumericType, D>(sphere, mesh, true, true).apply();

  auto writer = lsVTKWriter<NumericType>();
  writer.setMesh(mesh);
  writer.setFileName("curvatures.vtk");
  writer.apply();

  std::cout << "Finished" << std::endl;

  return 0;
}
