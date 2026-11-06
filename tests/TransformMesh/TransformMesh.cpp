#include <cmath>
#include <cstdio>

#include <lsMesh.hpp>
#include <lsTransformMesh.hpp>

/**
  TransformMesh with short vectors: a translation by (0, 0, -2e-4), as
  MakeHole applies to the cylinder that cuts its mask (twice its eps), a
  rotation about an axis given by a short vector and a scaling by a small
  factor move every node as given; a translation by the zero vector leaves
  the mesh as it is.
  \example TransformMesh.cpp
*/

namespace ls = viennals;

namespace {
ls::SmartPointer<ls::Mesh<double>> makeMesh() {
  auto mesh = ls::SmartPointer<ls::Mesh<double>>::New();
  mesh->insertNextNode({0., 0., 0.});
  mesh->insertNextNode({1., 0., 0.});
  mesh->insertNextNode({0., 1., 0.5});
  mesh->insertNextTriangle({0, 1, 2});
  return mesh;
}

int failures = 0;

void check(const char *what, const ls::SmartPointer<ls::Mesh<double>> &mesh,
           const double expected[3][3]) {
  for (int n = 0; n < 3; ++n)
    for (int i = 0; i < 3; ++i)
      if (std::abs(mesh->nodes[n][i] - expected[n][i]) > 1e-12) {
        std::printf("%s: node %d coordinate %d is %.15g, expected %.15g\n",
                    what, n, i, mesh->nodes[n][i], expected[n][i]);
        ++failures;
      }
}
} // namespace

int main() {
  {
    auto mesh = makeMesh();
    ls::TransformMesh<double>(mesh, ls::TransformEnum::TRANSLATION,
                              ls::Vec3D<double>{0., 0., -2e-4})
        .apply();
    const double expected[3][3] = {
        {0., 0., -2e-4}, {1., 0., -2e-4}, {0., 1., 0.5 - 2e-4}};
    check("translation by (0, 0, -2e-4)", mesh, expected);
  }
  {
    auto mesh = makeMesh();
    ls::TransformMesh<double>(mesh, ls::TransformEnum::TRANSLATION,
                              ls::Vec3D<double>{0., 0., 0.})
        .apply();
    const double expected[3][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.5}};
    check("translation by the zero vector", mesh, expected);
  }
  {
    // a quarter turn about the z axis, given by a vector of length 1e-4
    auto mesh = makeMesh();
    ls::TransformMesh<double>(mesh, ls::TransformEnum::ROTATION,
                              ls::Vec3D<double>{0., 0., 1e-4}, M_PI_2)
        .apply();
    const double expected[3][3] = {{0., 0., 0.}, {0., 1., 0.}, {-1., 0., 0.5}};
    check("rotation about (0, 0, 1e-4)", mesh, expected);
  }
  {
    auto mesh = makeMesh();
    ls::TransformMesh<double>(mesh, ls::TransformEnum::SCALE,
                              ls::Vec3D<double>{1e-4, 1e-4, 1e-4})
        .apply();
    const double expected[3][3] = {
        {0., 0., 0.}, {1e-4, 0., 0.}, {0., 1e-4, 0.5e-4}};
    check("scaling by 1e-4", mesh, expected);
  }
  if (failures > 0) {
    std::printf("TransformMesh: %d failures\n", failures);
    return 1;
  }
  std::printf("TransformMesh: all transforms applied\n");
  return 0;
}
