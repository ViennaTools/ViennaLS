#include <hrleDenseIterator.hpp>
#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write and read different
  meshes created by the algorithms lsToVoxelMesh and lsToSurfaceMesh.
  \example MakeSphere.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const double radius = 27.3;
  const ls::VectorType<double, D> centre(5., 0.);

  ls::MakeGeometry<double, 2>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
      .apply();

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;

  ls::ToMesh<double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "initial.vtk").apply();

  ls::Prune<double, D>(levelSet).apply();
  std::cout << "After prune: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  ls::ToMesh<double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "after_prune.vtk").apply();

  ls::Expand<double, D>(levelSet, 4).apply();
  std::cout << "After Expand: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  ls::ToMesh<double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "after_expand.vtk").apply();

  ls::Reduce<double, D>(levelSet, 2).apply();
  std::cout << "After Reduce: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "Sphere2D.vtk").apply();

  ls::ToMesh<double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "after_reduce.vtk").apply();

  ls::ToVoxelMesh<double, D>(levelSet, mesh).apply();
  mesh->print();
  ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTU, "Sphere.vtu").apply();

  return 0;
}
