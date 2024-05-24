#include <algorithm>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to use void detection.
  \example VoidDetection.cpp
*/

namespace ls = viennals;

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1.;
    return velocity;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<double, 3>({});
  }
};

int main() {
  constexpr int D = 2;
  omp_set_num_threads(1);

  double extent = 10;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {0., 1.};

  ls::MakeGeometry<double, D>(
      substrate, ls::SmartPointer<ls::Plane<double, D>>::New(origin, normal))
      .apply();
  {
    auto hole = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    origin[1] = -5.;
    ls::MakeGeometry<double, D>(
        hole, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, 3.))
        .apply();

    ls::BooleanOperation<double, D>(
        substrate, hole, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    auto explMesh = ls::SmartPointer<ls::Mesh<>>::New();

    std::cout << "Extracting..." << std::endl;
    ls::ToSurfaceMesh<double, D>(substrate, explMesh).apply();

    ls::VTKWriter<double>(explMesh, "before.vtp").apply();
  }

  ls::MarkVoidPoints<double, D>(substrate).apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "after.vtp").apply();
  }

  // Advection
  auto velocities = ls::SmartPointer<velocityField>::New();
  ls::Advect<double, D> advectionKernel(substrate, velocities);
  advectionKernel.setIgnoreVoids(true);
  advectionKernel.setSaveAdvectionVelocities(true);
  for (unsigned i = 0; i < 30; ++i) {
    {
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
      ls::VTKWriter<double>(mesh, "out-" + std::to_string(i) + ".vtp").apply();

      ls::MarkVoidPoints<double, D>(substrate).apply();
      ls::ToMesh<double, D>(substrate, mesh).apply();

      ls::VTKWriter<double>(mesh, "ls-out-" + std::to_string(i) + ".vtp")
          .apply();
    }
    advectionKernel.apply();
  }

  return 0;
}
