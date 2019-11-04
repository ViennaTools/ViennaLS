#include <iostream>

#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {

  constexpr int D = 2;

  omp_set_num_threads(1);

  // double extent = 100;
  // double gridDelta = 1;
  //
  // double bounds[2 * D] = {-extent, extent, -extent, extent};
  // lsDomain<double, D>::BoundaryType boundaryCons[D];
  // for (unsigned i = 0; i < D; ++i)
  //   boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  // lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);

  lsPointCloud<double, D> cloud;
  cloud.insertNextPoint(hrleVectorType<double, D>(-1, 0));
  cloud.insertNextPoint(hrleVectorType<double, D>(1, 0));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, -1));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, 1));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, 0));
  cloud.insertNextPoint(hrleVectorType<double, D>(0.5, 0.2));
  cloud.insertNextPoint(hrleVectorType<double, D>(1.2, 0.8));
  cloud.insertNextPoint(hrleVectorType<double, D>(0.1, 0.5));
  cloud.insertNextPoint(hrleVectorType<double, D>(-1, 0.2));

  lsMesh pointMesh;
  for (unsigned i = 0; i < cloud.points.size(); ++i) {
    pointMesh.nodes.push_back(
        hrleVectorType<double, 3>(cloud.points[i][0], cloud.points[i][1], 0.));
    pointMesh.vertices.push_back(hrleVectorType<unsigned, 1>(i));
  }
  lsVTKWriter(pointMesh).writeVTKLegacy("points.vtk");

  lsMesh mesh;
  lsConvexHull<double, D>(mesh, cloud).apply();

  mesh.print();
  lsVTKWriter(mesh).writeVTKLegacy("hull.vtk");

  // {
  //   lsMesh mesh;
  //   std::cout << "Extracting..." << std::endl;
  //   lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   mesh.print();
  //   lsVTKWriter(mesh).writeVTKLegacy("after2D.vtk");
  // }

  return 0;
}
