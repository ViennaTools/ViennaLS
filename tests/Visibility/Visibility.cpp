#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <lsCalculateVisibilities.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {

  constexpr int D = 2;

  double gridDelta = 0.4;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  ls::Expand<double, D>(sphere1, 5).apply();
  ls::CalculateVisibilities<double, D>(sphere1, ls::Vec3D<double>{0., -1.0, 0.})
      .apply();

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(sphere1, mesh).apply();
  ls::VTKWriter<double>(mesh, "visibility_test.vtp").apply();

  // a Euclidean level set holds distances far from the surface: a node
  // outside is shadowed exactly where its ray towards the source passes
  // through the solid, below the circle, and seen elsewhere, also beside the
  // circle, where the ray passes close to it
  auto circle = ls::SmartPointer<ls::Domain<double, D>>::New(gridDelta);
  circle->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  circle->setLevelSetWidth(12);
  ls::MakeGeometry<double, D>(
      circle, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  ls::CalculateVisibilities<double, D>(circle, ls::Vec3D<double>{0., -1.0, 0.})
      .apply();
  const auto *visibilities =
      circle->getPointData().getScalarData("Visibilities");
  int checked = 0, wrong = 0;
  for (viennahrle::ConstSparseIterator<ls::Domain<double, D>::DomainType> it(
           circle->getDomain());
       !it.isFinished(); it.next()) {
    if (!it.isDefined() || it.getValue() < 0.)
      continue;
    const double x = it.getStartIndices(0) * gridDelta - origin[0];
    const double y = it.getStartIndices(1) * gridDelta - origin[1];
    // rays that graze the circle
    if (std::abs(std::abs(x) - radius) < 2. * gridDelta)
      continue;
    const bool shadowed = std::abs(x) < radius && y < 0.;
    ++checked;
    if ((*visibilities)[it.getPointId()] != (shadowed ? 0. : 1.))
      ++wrong;
  }
  std::cout << "Euclidean circle: " << wrong << " of " << checked
            << " nodes outside with the wrong visibility" << std::endl;
  if (checked == 0 || wrong > 0)
    return EXIT_FAILURE;

  return 0;
}