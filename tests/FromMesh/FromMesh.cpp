#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include <lsCompareVolume.hpp>
#include <lsFromMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

// Signed distances of a sphere in a band wider than the 2.5 grid spacings
// FromMesh keeps, in random order, as exported by other level set tools.
template <int D> void testWideUnsortedBand() {
  double extent = 15;
  double gridDelta = 0.5;
  double radius = 7.3;

  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    bounds[2 * i] = -extent;
    bounds[2 * i + 1] = extent;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }

  std::vector<std::pair<std::array<double, 3>, double>> band;
  const int n = static_cast<int>(extent / gridDelta);
  const int nz = (D == 3) ? n : 0;
  for (int k = -nz; k <= nz; ++k) {
    for (int j = -n; j <= n; ++j) {
      for (int i = -n; i <= n; ++i) {
        std::array<double, 3> x{i * gridDelta, j * gridDelta, k * gridDelta};
        double value =
            (std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - radius) /
            gridDelta;
        if (std::abs(value) <= 3.)
          band.emplace_back(x, value);
      }
    }
  }
  std::shuffle(band.begin(), band.end(), std::mt19937(42));

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  std::vector<double> values;
  for (const auto &[x, value] : band) {
    mesh->insertNextNode(x);
    values.push_back(value);
  }
  mesh->getPointData().insertNextScalarData(values, "LSValues");

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  ls::FromMesh<double, D>(levelSet, mesh).apply();
  LSTEST_ASSERT_VALID_LS(levelSet, double, D)

  auto reference = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  double origin[D] = {};
  ls::MakeGeometry<double, D>(
      reference, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  ls::CompareVolume<double, D> compare(reference, levelSet);
  compare.apply();
  VC_TEST_ASSERT(compare.getVolumeMismatch() < 1e-6)
}

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  std::cout << "Writing" << std::endl;
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere.vtk").apply();
  }

  std::cout << "Reading" << std::endl;
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::VTKReader<double>(mesh, "sphere.vtk").apply();
    auto newLS = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    ls::FromMesh<double, D>(newLS, mesh).apply();

    std::cout << "Writing new" << std::endl;
    auto newMesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(newLS, newMesh).apply();
    ls::VTKWriter<double>(mesh, "newMesh.vtk").apply();
  }

  testWideUnsortedBand<2>();
  testWideUnsortedBand<3>();

  return 0;
}
