#include <iostream>
#include <algorithm>


#include <lsGraph.hpp>

/**
  Example showing how to use void detection.
  \example AdvectionPlane.cpp
*/

int main() {

  lsInternal::lsGraph graph;

  graph.insertNextEdge(0, 1);

  for(unsigned i = 0; i < 100; ++i) {
    if(i%10) {
      graph.insertNextEdge(i, i+1);
    }
  }

  auto result = graph.getConnectedComponents();
  auto it = std::unique(result.begin(), result.end());
  result.resize(std::distance(result.begin(), it));
  for(auto it = result.begin(); it != result.end(); ++it){
    std::cout << *it << std::endl;
  }


  // constexpr int D = 2;
  // omp_set_num_threads(1);
  //
  // double extent = 25;
  // double gridDelta = 1;
  //
  // double bounds[2 * D] = {-extent, extent, -extent, extent};
  // lsDomain<double, D>::BoundaryType boundaryCons[D];
  // for (unsigned i = 0; i < D; ++i)
  //   boundaryCons[i] = lsDomain<double, D>::BoundaryType::SYMMETRIC_BOUNDARY;
  // lsDomain<double, D> plane(bounds, boundaryCons, gridDelta);
  //
  // double origin[D] = {0., 0.};
  // double normal[D] = {1., 1.};
  //
  // lsMakeGeometry<double, D>(plane).makePlane(origin, normal);
  // {
  //   lsMesh mesh;
  //   lsMesh explMesh;
  //
  //   std::cout << "Extracting..." << std::endl;
  //   lsToExplicitMesh<double, D>(plane, explMesh).apply();
  //   lsToMesh<double, D>(plane, mesh).apply();
  //
  //   mesh.print();
  //   lsVTKWriter(explMesh).writeVTKLegacy("before.vtk");
  //   lsVTKWriter(mesh).writeVTKLegacy("beforeLS.vtk");
  // }
  //
  // // fill vector with lsDomain pointers
  // std::vector<lsDomain<double, D> *> lsDomains;
  // lsDomains.push_back(&plane);
  //
  // velocityField velocities;
  //
  // std::cout << "number of Points: " << plane.getDomain().getNumberOfPoints()
  //           << std::endl;
  //
  // std::cout << "Advecting" << std::endl;
  // lsAdvect<double, D> advectionKernel(lsDomains, velocities);
  // advectionKernel.apply();
  // double advectionTime = advectionKernel.getAdvectionTime();
  // std::cout << "Time difference: " << advectionTime << std::endl;
  //
  // lsPrune<double, D>(plane).apply();
  // lsExpand<double, D>(plane, 2).apply();
  //
  // std::cout << "Extracting..." << std::endl;
  // lsMesh mesh;
  // lsToExplicitMesh<double, D>(plane, mesh).apply();
  //
  // // mesh.print();
  //
  // lsVTKWriter(mesh).writeVTKLegacy("after.vtk");

  return 0;
}
