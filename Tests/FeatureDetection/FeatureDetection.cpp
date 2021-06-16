#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsDetectFeatures.hpp>

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

  NumericType gridDelta = 4;

  std::vector<NumericType> planeNormal{0, 0, 1};
  std::vector<NumericType> origin{0, 0, 0.3};

  std::cout << "Creating trench..." << std::endl;

  double extent = 47.3;
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
                                  (NumericType)-extent - 1, -49.8};
      NumericType maxCorner[3] = {(NumericType)(extent / 4.),
                                  (NumericType)extent + 1, 5.0};
      auto box =
          lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
      lsMakeGeometry<NumericType, D>(trench, box).apply();
    } else {
      NumericType minCorner[2] = {(NumericType)(-extent / 4.), -49.8};
      NumericType maxCorner[2] = {(NumericType)(extent / 4.), 5.0};
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

  lsDetectFeatures<NumericType, D> CurvatureFlagger(
      levelSet, 1e-3, lsFeatureDetectionEnum::CURVATURE, "Features_Curve");

  CurvatureFlagger.apply();

  std::cout << "Flagging Normals..." << std::endl;

  lsDetectFeatures<NumericType, D> NormalsFlagger(
      levelSet, 1e-3, lsFeatureDetectionEnum::NORMALS_ANGLE, "Features_Angle");

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