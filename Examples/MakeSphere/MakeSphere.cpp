#include <hrleIterator.hpp>
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
#include <lsVTKWriter.hpp>

void printLS(lsDomain_double_2 &domain) {
  hrleIterator<hrleDomain<double, 2>> it(domain.getDomain());
  int y = domain.getDomain().getMinRunBreak(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      std::cout << std::endl << std::endl;
    }
    if (it.getValue() == lsDomain_double_2::POS_VALUE) {
      std::cout << std::setw(8) << "+oo"
                << " ";
    } else if (it.getValue() == lsDomain_double_2::NEG_VALUE) {
      std::cout << std::setw(8) << "-oo"
                << " ";
    } else {
      std::cout << std::setw(8) << std::setprecision(2) << it.getValue() << " ";
    }
    ++it;
  }
  std::cout << std::endl;
}

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  lsDomain_double_2 levelSet;

  const double radius = 27.3;
  const hrleVectorType<int, D> centre(5., 0.); // all zeros

  lsMakeGeometry<double, 2>(levelSet).makeSphere(centre, radius);

  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;

  printLS(levelSet);

  lsPrune<double, D>(levelSet).apply();
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  printLS(levelSet);

  lsExpand<double, D>(levelSet).apply(4);
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  printLS(levelSet);

  lsReduce<double, D>(levelSet).apply(2);
  std::cout << "Number of points: " << levelSet.getDomain().getNumberOfPoints()
            << std::endl;
  std::cout << "Width: " << levelSet.getLevelSetWidth() << std::endl;

  // printLS(levelSet);
  lsMesh mesh;

  lsToExplicitMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("Sphere2D.vtk");

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh).writeVTKLegacy("SphereLS.vtk");

  return 0;
}
