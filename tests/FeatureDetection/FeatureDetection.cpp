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

namespace ls = viennals;

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
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  {
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(levelSet, plane).apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;

    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    if constexpr (D == 3) {
      NumericType minCorner[3] = {(NumericType)(-extent / 4.),
                                  (NumericType)-extent - 1, -49.8};
      NumericType maxCorner[3] = {(NumericType)(extent / 4.),
                                  (NumericType)extent + 1, 5.0};
      auto box =
          ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
      ls::MakeGeometry<NumericType, D>(trench, box).apply();
    } else {
      NumericType minCorner[2] = {(NumericType)(-extent / 4.), -49.8};
      NumericType maxCorner[2] = {(NumericType)(extent / 4.), 5.0};
      auto box =
          ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
      ls::MakeGeometry<NumericType, D>(trench, box).apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    ls::BooleanOperation<NumericType, D>(
        levelSet, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  levelSet->getDomain().segment();

  std::cout << "Expanding..." << std::endl;

  ls::Expand<NumericType, D> expander(levelSet, 5);

  expander.apply();

  std::cout << "Flagging Curvatures..." << std::endl;

  ls::DetectFeatures<NumericType, D> CurvatureFlagger(
      levelSet, 1e-3, ls::FeatureDetectionEnum::CURVATURE);

  CurvatureFlagger.apply();

  // now rename markers, so they are not overwritten by next calls
  auto &pointData = levelSet->getPointData();
  for (unsigned i = 0; i < pointData.getScalarDataSize(); ++i) {
    if (pointData.getScalarDataLabel(i) == "FeatureMarkers") {
      pointData.setScalarDataLabel(i, "FeatureMarkers_Curvature");
    }
  }

  std::cout << "Flagging Normals..." << std::endl;

  ls::DetectFeatures<NumericType, D> NormalsFlagger(
      levelSet, 1e-3, ls::FeatureDetectionEnum::NORMALS_ANGLE);

  NormalsFlagger.apply();

  std::cout << "Writing Output..." << std::endl;

  auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  ls::ToMesh<NumericType, D>(levelSet, mesh, true, true).apply();

  auto writer = ls::VTKWriter<NumericType>();
  writer.setMesh(mesh);
  writer.setFileName("Features.vtk");
  writer.apply();

  std::cout << "Finished" << std::endl;

  return 0;
}