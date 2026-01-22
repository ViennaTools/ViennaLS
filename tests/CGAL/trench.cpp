#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDelaunay2D.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

using NumericType = double;

// implement own velocity field
class velocityField : public ls::VelocityField<NumericType> {
public:
  NumericType
  getScalarVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    NumericType velocity = 1.;
    return velocity;
  }
};

int main() {

  constexpr int D = 2;

  NumericType extent = 30;
  NumericType gridDelta = 0.5;

  double bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -extent;
    bounds[2 * i + 1] = extent;
  }

  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  NumericType origin[D] = {0.};
  NumericType planeNormal[D] = {0.};
  planeNormal[D - 1] = 1.;

  {
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    NumericType ylimit = extent / 3.;
    NumericType minCorner[D];
    NumericType maxCorner[D];
    if constexpr (D == 2) {
      minCorner[0] = -ylimit;
      maxCorner[0] = ylimit;
    } else {
      minCorner[0] = -extent - 1;
      maxCorner[0] = extent + 1;
      minCorner[1] = -ylimit;
      maxCorner[1] = ylimit;
    }
    minCorner[D - 1] = -15.;
    maxCorner[D - 1] = 1.;
    auto box =
        ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<NumericType, D>(trench, box).apply();
    ls::BooleanOperation<NumericType, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  auto newLayer = ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);
  auto velocities = ls::SmartPointer<velocityField>::New();
  ls::Advect<NumericType, D> advectionKernel;
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setAdvectionTime(4.);
  advectionKernel.apply();

  auto mesh = ls::Mesh<NumericType>::New();

  ls::Delaunay2D<NumericType> delaunay;
  delaunay.setMesh(mesh);
  delaunay.insertNextLevelSet(substrate);
  delaunay.insertNextLevelSet(newLayer);
  delaunay.setMaxTriangleSize(gridDelta * 2.);
  delaunay.setBottomExtent(10);
  delaunay.apply();

  ls::VTKWriter<NumericType>(mesh, "trench").apply();

  return 0;
}
