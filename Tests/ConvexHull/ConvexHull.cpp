#include <iostream>
#include <random>

#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#define PI 3.141592

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

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0., 1.);

  auto cloud = lsSmartPointer<lsPointCloud<double, D>>::New();

  // generate a circle of points
  int numberOfPoints = 10000;
  double radius = 2.;
  for (int i = 0; i < numberOfPoints; ++i) {
    double angle = dis(gen) * 2 * PI;
    double r = radius * sqrt(dis(gen));

    double x = r * cos(angle);
    double y = r * sin(angle);
    cloud->insertNextPoint(hrleVectorType<double, D>(x, y));
  }
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, -1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, 1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0.5, 0.2));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1.2, 0.8));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0.1, 0.5));
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, 0.2));

  auto pointMesh = lsSmartPointer<lsMesh<>>::New();
  for (unsigned i = 0; i < cloud->points.size(); ++i) {
    pointMesh->nodes.push_back(
        std::array<double, 3>{cloud->points[i][0], cloud->points[i][1], 0.});
    pointMesh->vertices.push_back(std::array<unsigned, 1>{i});
  }
  lsVTKWriter<double>(pointMesh, lsFileFormatEnum::VTP, "points.vtp").apply();

  lsMakeGeometry<double, D> geom;
  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  geom.setLevelSet(levelSet);
  geom.setGeometry(cloud);
  geom.apply();

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
    mesh->print();
    lsVTKWriter<double>(mesh, "LSMesh.vtk").apply();
  }

  return 0;
}
