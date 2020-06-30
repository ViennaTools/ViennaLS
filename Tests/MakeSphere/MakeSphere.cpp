#include <hrleDenseIterator.hpp>
#include <hrleVectorType.hpp>
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

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh>::New();

  const double radius = 27.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, 2>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "initial.vtk").apply();

  lsPrune<double, D>(levelSet).apply();
  std::cout << "After prune: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "after_prune.vtk").apply();

  lsExpand<double, D>(levelSet, 4).apply();
  std::cout << "After Expand: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "after_expand.vtk").apply();

  lsReduce<double, D>(levelSet, 2).apply();
  std::cout << "After Reduce: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet->getLevelSetWidth() << std::endl;

  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "Sphere2D.vtk").apply();

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "after_reduce.vtk").apply();

  lsToVoxelMesh<double, D>(levelSet, mesh).apply();
  mesh->print();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "Sphere.vtu").apply();

  return 0;
}
