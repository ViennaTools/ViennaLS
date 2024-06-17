#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

using NumericType = float;
constexpr int D = 2;

namespace ls = viennals;

void outputDomain(ls::SmartPointer<ls::Domain<NumericType, D>> domain,
                  std::string fileName) {
  auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  ls::ToSurfaceMesh(domain, mesh).apply();
  ls::VTKWriter(mesh, ls::FileFormatEnum::VTP, fileName + "_surface.vtp")
      .apply();

  ls::ToMesh(domain, mesh).apply();
  ls::VTKWriter(mesh, ls::FileFormatEnum::VTP, fileName + "_points.vtp")
      .apply();
}

void makeGeometry(ls::SmartPointer<ls::Domain<NumericType, D>> domain) {
  auto plane =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(domain->getGrid());

  NumericType origin[D] = {0., 0.35};
  NumericType normal[D] = {0., 1.};

  ls::MakeGeometry(
      plane, ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, normal))
      .apply();

  ls::BooleanOperation(domain, plane, ls::BooleanOperationEnum::UNION).apply();

  origin[0] = -8.;
  origin[1] = -9.5;
  double radius = 5.1;
  for (unsigned i = 0; i < 2; ++i) {
    auto hole =
        ls::SmartPointer<ls::Domain<NumericType, D>>::New(domain->getGrid());

    ls::MakeGeometry(
        hole, ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
        .apply();

    ls::BooleanOperation(domain, hole,
                         ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
    origin[0] += 16.;
  }
}

int main() {
  omp_set_num_threads(1);

  NumericType extent = 15;
  NumericType gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::BoundaryConditionEnum<D> boundaryCons[D];
  boundaryCons[0] = ls::BoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::BoundaryConditionEnum<D>::INFINITE_BOUNDARY;

  auto domain = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  makeGeometry(domain);

  // outputDomain(domain, "initial");

  auto marker = ls::MarkVoidPoints(domain);
  marker.setSaveComponentIds(true);
  marker.apply();

  auto markers = domain->getPointData().getScalarData("VoidPointMarkers");
  LSTEST_ASSERT(markers)

  // check if void points are set correctly
  bool correct = true;
  for (auto it = hrleConstSparseIterator<
           typename ls::Domain<NumericType, D>::DomainType>(
           domain->getDomain());
       !it.isFinished(); ++it) {
    // skip undefined runs
    if (!it.isDefined())
      continue;

    if (it.getStartIndices()[1] < -2) {
      if (auto m = markers->at(it.getPointId()); m == 0) {
        std::cout << "ERROR: Wrong VoidPointMarker " << m << " at "
                  << it.getStartIndices() << std::endl;
        correct = false;
      }
    } else {
      if (auto m = markers->at(it.getPointId()); m != 0) {
        std::cout << "ERROR: Wrong VoidPointMarker " << m << " at "
                  << it.getStartIndices() << std::endl;
        correct = false;
      }
    }
  }
  LSTEST_ASSERT(correct)

  // std::cout << "voidMarkers: " <<
  // domain->getPointData().getScalarData("VoidPointMarkers")->size() <<
  // std::endl; std::cout << "componentMarkers: " <<
  // domain->getPointData().getScalarData("ConnectedComponentId")->size() <<
  // std::endl;

  // outputDomain(domain, "marked");

  return 0;
}
