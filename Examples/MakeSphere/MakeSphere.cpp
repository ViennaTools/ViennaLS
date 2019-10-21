#include <hrleDenseIterator.hpp>
#include <hrleVectorType.hpp>
#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write and read different
  meshes created by the algorithms lsToVoxelMesh and lsToExplicitMesh.
  \example MakeSphere.cpp
*/

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  lsDomain<double, D> levelSet;
  lsMesh mesh;

  const double radius = 27.3;
  const hrleVectorType<int, D> centre(5., 0.); // all zeros

  lsMakeGeometry<double, 2>(levelSet).makeSphere(centre, radius);

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("initial.vtk");

  lsPrune<double, D>(levelSet).apply();
  std::cout << "After prune: " << std::endl;
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("after_prune.vtk");

  lsExpand<double, D>(levelSet).apply(4);
  std::cout << "After Expand: " << std::endl;
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("after_expand.vtk");

  lsReduce<double, D>(levelSet).apply(2);
  std::cout << "After Reduce: " << std::endl;
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  lsToExplicitMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("Sphere2D.vtk");

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("after_reduce.vtk");

  lsToVoxelMesh<double, D>(levelSet, mesh).apply();
  mesh.print();
  lsVTKWriter(mesh).writeVTU("Sphere.vtu");

  return 0;
}
