#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsFeatureDetection.hpp>

#include <omp.h>

/**
  Minimal example of how to use feature detection on a level set function.
  Works in 2D and 3D.
  Outputs the detected features into a .vtk file.
  \example calculateCurvatures.cpp
*/

constexpr int D = 3;
typedef double NumericType;

int main() {

  omp_set_num_threads(4);

  NumericType gridDelta = 0.5;

  std::vector<NumericType> planeNormal;
  std::vector<NumericType> origin;

  for (int i = 0; i < D; i++) {
    planeNormal.push_back(0.);
    origin.push_back(0.);
  }
  planeNormal[D - 1] = 1.;

  std::cout << "Creating trench..." << std::endl;

  double extent = 50;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  {
    auto plane =
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(levelSet, plane).apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;

    auto trench = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    if (D == 3) {
      NumericType minCorner[3] = {(NumericType)(-extent / 4.),
                                  (NumericType)-extent - 1, -50.};
      NumericType maxCorner[3] = {(NumericType)(extent / 4.),
                                  (NumericType)extent + 1, 1.0};
      auto box =
          lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
      lsMakeGeometry<NumericType, D>(trench, box).apply();
    } else {
      NumericType minCorner[2] = {(NumericType)(-extent / 4.), -50.};
      NumericType maxCorner[2] = {(NumericType)(extent / 4.), 1.0};
      auto box =
          lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
      lsMakeGeometry<NumericType, D>(trench, box).apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    lsBooleanOperation<NumericType, D>(
        levelSet, trench, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  levelSet->getDomain().segment();

  std::cout << "Expanding..." << std::endl;

  lsExpand<NumericType, D> expander(levelSet, 5);

  expander.apply();

  std::cout << "Flagging Curvatures..." << std::endl;

  lsFeatureDetection<NumericType, D> CurvatureFlagger(levelSet, 1e-3);

  CurvatureFlagger.apply();

  std::cout << "Flagging Normals..." << std::endl;

  lsFeatureDetection<NumericType, D> NormalsFlagger(
      levelSet, 0.16, FeatureDetectionMethod::NORMALS_ANGLE, "Features_Angle");

  NormalsFlagger.apply();

  std::cout << "Writing Output..." << std::endl;

  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  lsToMesh<NumericType, D>(levelSet, mesh, true, true).apply();

  auto writer = lsVTKWriter<NumericType>();
  writer.setMesh(mesh);
  writer.setFileName("Features.vtk");
  writer.apply();

  std::cout << "Finished" << std::endl;

  return 0;
}