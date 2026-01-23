#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsDelaunay3D.hpp>

namespace ls = viennals;

int main() {
  constexpr int D = 3;

  omp_set_num_threads(4);

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const double radius = 10.0;
  const ls::VectorType<double, D> centre{0., 0., 0.};

  ls::MakeGeometry<double, 3>(levelSet,
                              ls::Sphere<double, D>::New(centre, radius))
      .apply();

  ls::Delaunay3D<double> delaunay;
  delaunay.insertNextLevelSet(levelSet);
  delaunay.setMesh(mesh);
  delaunay.apply();

  return 0;
}
