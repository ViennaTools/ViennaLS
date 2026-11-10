#include <cmath>
#include <iostream>
#include <limits>

#include <lsCompareSurfaceDistance.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>

#include <vcTestAsserts.hpp>

/**
  Test for lsCompareSurfaceDistance, which measures two surfaces from their
  uniformly sampled elements to the nearest point of the other surface. The
  surfaces are meshes or level sets of any grid spacing and normalization.
  - Concentric circles (2D) and spheres (3D) as meshes with different vertices:
    every distance is the difference of the radii.
  - A sphere as Manhattan and Euclidean level sets on two grid spacings
    against the exact sphere, and the two grids against each other.
  - A box off the grid as a Euclidean level set, whose dual contouring surface
    is the exact box, and as a Manhattan level set, whose edges are cut.
  - The result does not change when the sampling spacing is halved.
  - A box restricts the samples, and a missing surface gives NaN.
  \example CompareSurfaceDistance.cpp
*/

namespace ls = viennals;

using T = double;

// a circle (2D) or a UV sphere (3D) of the given radius around the origin
template <int D>
ls::SmartPointer<ls::Mesh<T>> sphereMesh(T radius, int segments) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  if constexpr (D == 2) {
    for (int k = 0; k < segments; ++k) {
      const T a = 2. * M_PI * k / segments;
      mesh->insertNextNode({radius * std::cos(a), radius * std::sin(a), 0.});
      mesh->insertNextLine({unsigned(k), unsigned((k + 1) % segments)});
    }
  } else {
    const int rings = segments / 2;
    mesh->insertNextNode({0., 0., radius});
    for (int i = 1; i < rings; ++i) {
      const T theta = M_PI * i / rings;
      for (int j = 0; j < segments; ++j) {
        const T phi = 2. * M_PI * j / segments;
        mesh->insertNextNode({radius * std::sin(theta) * std::cos(phi),
                              radius * std::sin(theta) * std::sin(phi),
                              radius * std::cos(theta)});
      }
    }
    const unsigned south = mesh->insertNextNode({0., 0., -radius});
    auto ring = [&](int i, int j) {
      return unsigned(1 + (i - 1) * segments + (j % segments));
    };
    for (int j = 0; j < segments; ++j) {
      mesh->insertNextTriangle({0u, ring(1, j), ring(1, j + 1)});
      mesh->insertNextTriangle(
          {south, ring(rings - 1, j + 1), ring(rings - 1, j)});
    }
    for (int i = 1; i + 1 < rings; ++i)
      for (int j = 0; j < segments; ++j) {
        mesh->insertNextTriangle(
            {ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)});
        mesh->insertNextTriangle(
            {ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)});
      }
  }
  return mesh;
}

// the surface of an axis-aligned box
template <int D>
ls::SmartPointer<ls::Mesh<T>> boxMesh(const T *lo, const T *hi) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  if constexpr (D == 2) {
    mesh->insertNextNode({lo[0], lo[1], 0.});
    mesh->insertNextNode({hi[0], lo[1], 0.});
    mesh->insertNextNode({hi[0], hi[1], 0.});
    mesh->insertNextNode({lo[0], hi[1], 0.});
    for (unsigned k = 0; k < 4; ++k)
      mesh->insertNextLine({k, (k + 1) % 4});
  } else {
    for (unsigned c = 0; c < 8; ++c)
      mesh->insertNextNode({(c & 1) ? hi[0] : lo[0], (c & 2) ? hi[1] : lo[1],
                            (c & 4) ? hi[2] : lo[2]});
    const unsigned quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4},
                                  {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
    for (const auto &q : quads) {
      mesh->insertNextTriangle({q[0], q[1], q[2]});
      mesh->insertNextTriangle({q[0], q[2], q[3]});
    }
  }
  return mesh;
}

template <int D>
ls::SmartPointer<ls::Domain<T, D>> newDomain(T gridDelta, bool euclidean) {
  T bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -10.;
    bounds[2 * i + 1] = 10.;
    boundaryCons[i] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;
  }
  auto domain =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  if (euclidean) {
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
  }
  return domain;
}

template <int D>
ls::SmartPointer<ls::Domain<T, D>> sphereLevelSet(T gridDelta, bool euclidean,
                                                  T radius) {
  auto domain = newDomain<D>(gridDelta, euclidean);
  T origin[D] = {};
  ls::MakeGeometry<T, D>(
      domain, ls::SmartPointer<ls::Sphere<T, D>>::New(origin, radius))
      .apply();
  return domain;
}

template <int D> void runTest() {
  std::cout << "== " << D << "D" << std::endl;
  const T radius = 5.;

  // concentric surfaces with different vertices: every distance is 0.3
  {
    const int segments = D == 2 ? 2000 : 180;
    auto inner = sphereMesh<D>(radius, segments);
    auto outer = sphereMesh<D>(radius + 0.3, segments + (D == 2 ? 300 : 20));
    ls::CompareSurfaceDistance<T, D> compare(inner, outer);
    if constexpr (D == 3)
      compare.setSamplingSpacing(0.2);
    compare.apply();
    // the facets of the meshes lie within this of their circles or spheres
    const T facet = D == 2 ? 1e-5 : 1.5e-3;
    std::cout << "concentric: forward " << compare.getForwardDistance()
              << ", backward " << compare.getBackwardDistance() << ", Chamfer "
              << compare.getChamferDistance() << ", RMS "
              << compare.getRMSChamferDistance() << ", maximum "
              << compare.getMaxDistance() << ", samples "
              << compare.getNumTargetPoints() << " / "
              << compare.getNumSamplePoints() << std::endl;
    VC_TEST_ASSERT(std::abs(compare.getForwardDistance() - 0.3) < 2 * facet);
    VC_TEST_ASSERT(std::abs(compare.getBackwardDistance() - 0.3) < 2 * facet);
    VC_TEST_ASSERT(std::abs(compare.getChamferDistance() - 0.3) < 2 * facet);
    VC_TEST_ASSERT(std::abs(compare.getRMSChamferDistance() - 0.3) < 2 * facet);
    VC_TEST_ASSERT(std::abs(compare.getMaxDistance() - 0.3) < 2 * facet);
  }

  // a sphere as level sets of both normalizations on two grid spacings
  auto exact = sphereMesh<D>(radius, D == 2 ? 4000 : 360);
  const T exactFacet = D == 2 ? 2e-6 : 4e-4;
  for (bool euclidean : {false, true})
    for (T gridDelta : {0.5, 0.25}) {
      auto sphere = sphereLevelSet<D>(gridDelta, euclidean, radius);
      ls::CompareSurfaceDistance<T, D> compare;
      compare.setLevelSetSample(sphere);
      compare.setMeshTarget(exact);
      compare.apply();
      std::cout << (euclidean ? "Euclidean" : "Manhattan") << " sphere at "
                << gridDelta << ": Chamfer " << compare.getChamferDistance()
                << ", maximum " << compare.getMaxDistance() << std::endl;
      VC_TEST_ASSERT(std::isfinite(compare.getChamferDistance()));
      VC_TEST_ASSERT(compare.getChamferDistance() < 0.1 * gridDelta);
      VC_TEST_ASSERT(compare.getMaxDistance() < 0.5 * gridDelta + exactFacet);
    }
  {
    // two grid spacings against each other
    auto coarse = sphereLevelSet<D>(0.5, false, radius);
    auto fine = sphereLevelSet<D>(0.25, false, radius);
    ls::CompareSurfaceDistance<T, D> compare(coarse, fine);
    compare.apply();
    std::cout << "Manhattan sphere, 0.5 against 0.25: Chamfer "
              << compare.getChamferDistance() << std::endl;
    VC_TEST_ASSERT(std::isfinite(compare.getChamferDistance()) &&
                   compare.getChamferDistance() < 0.05);
  }

  // a box off the grid: exact as a Euclidean level set, cut at the edges as a
  // Manhattan one
  const T lo[3] = {-3.3, -1.9, -2.6}, hi[3] = {2.7, 4.1, 1.8};
  auto boxExact = boxMesh<D>(lo, hi);
  T manhattanChamfer = 0.;
  for (bool euclidean : {true, false}) {
    auto box = newDomain<D>(0.5, euclidean);
    ls::MakeGeometry<T, D>(box, ls::SmartPointer<ls::Box<T, D>>::New(lo, hi))
        .apply();
    ls::CompareSurfaceDistance<T, D> compare;
    compare.setMeshTarget(boxExact);
    compare.setLevelSetSample(box);
    auto output = ls::SmartPointer<ls::Mesh<T>>::New();
    compare.setOutputMeshSample(output);
    compare.apply();
    std::cout << (euclidean ? "Euclidean" : "Manhattan") << " box: Chamfer "
              << compare.getChamferDistance() << ", maximum "
              << compare.getMaxDistance() << std::endl;
    if (euclidean) {
      VC_TEST_ASSERT(compare.getChamferDistance() < 1e-9);
      VC_TEST_ASSERT(compare.getMaxDistance() < 1e-9);
    } else {
      manhattanChamfer = compare.getChamferDistance();
      VC_TEST_ASSERT(manhattanChamfer > 1e-3);
    }
    VC_TEST_ASSERT(output->getNodes().size() > 0 &&
                   output->getPointData().getScalarData("DistanceToTarget") !=
                       nullptr);
  }

  // halving the sampling spacing changes the result by less than 1 %
  {
    auto box = newDomain<D>(0.5, false);
    ls::MakeGeometry<T, D>(box, ls::SmartPointer<ls::Box<T, D>>::New(lo, hi))
        .apply();
    T chamfer[2];
    for (int k = 0; k < 2; ++k) {
      ls::CompareSurfaceDistance<T, D> compare;
      compare.setMeshTarget(boxExact);
      compare.setLevelSetSample(box);
      compare.setSamplingSpacing(k == 0 ? 0.1 : 0.05);
      compare.apply();
      chamfer[k] = compare.getChamferDistance();
    }
    std::cout << "sampling 0.1 and 0.05: Chamfer " << chamfer[0] << ", "
              << chamfer[1] << std::endl;
    VC_TEST_ASSERT(std::abs(chamfer[0] - chamfer[1]) < 0.01 * chamfer[1]);
    VC_TEST_ASSERT(std::abs(chamfer[1] - manhattanChamfer) <
                   0.01 * manhattanChamfer);
  }

  // a box around one face restricts the samples
  {
    auto inner = sphereMesh<D>(radius, D == 2 ? 400 : 60);
    auto outer = sphereMesh<D>(radius + 0.3, D == 2 ? 500 : 80);
    ls::CompareSurfaceDistance<T, D> all(inner, outer), half(inner, outer);
    std::array<T, D> boxMin, boxMax;
    for (int i = 0; i < D; ++i) {
      boxMin[i] = -10.;
      boxMax[i] = 10.;
    }
    boxMin[D - 1] = 0.;
    half.setBox(boxMin, boxMax);
    all.apply();
    half.apply();
    std::cout << "box: samples " << half.getNumTargetPoints() << " of "
              << all.getNumTargetPoints() << std::endl;
    VC_TEST_ASSERT(half.getNumTargetPoints() > 0 &&
                   half.getNumTargetPoints() < all.getNumTargetPoints());
    VC_TEST_ASSERT(std::abs(half.getChamferDistance() - 0.3) < 0.01);
  }

  // a missing surface
  {
    ls::CompareSurfaceDistance<T, D> compare;
    compare.setMeshTarget(exact);
    compare.apply();
    VC_TEST_ASSERT(std::isnan(compare.getChamferDistance()));
  }
}

int main() {
  omp_set_num_threads(4);
  runTest<2>();
  runTest<3>();
  return 0;
}
