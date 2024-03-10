#include <iostream>
#include <random>
// #include <cmath>
// #include <cstdlib>
// #include <ctime>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#define PI 3.141592

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  // double extent = 100;
  // double gridDelta = 1;
  //
  // double bounds[2 * D] = {-extent, extent, -extent, extent};
  // lsDomain<double, D>::BoundaryType boundaryCons[D];
  // for (unsigned i = 0; i < D; ++i)
  //   boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  // lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);

  // auto seed = 1572962798; // std::time(nullptr);
  // std::cout << seed << std::endl;
  // std::srand(seed);

  auto cloud = lsSmartPointer<lsPointCloud<double, D>>::New();

  // generate a circle of points ----------------------------------------------
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> dis(0., 1.);
  // int numberOfPoints = 10000;
  // double radius = 4.;
  // for (int i = 0; i < numberOfPoints; ++i) {
  //   // double theta = std::rand() * PI;
  //   // double phi = std::rand() * 2 * PI;
  //   // double r = radius * sqrt(std::rand());
  //
  //   double theta = dis(gen) * PI;
  //   double phi = dis(gen) * 2 * PI;
  //   double r = radius * sqrt(dis(gen));
  //
  //   double x = r * sin(theta) * cos(phi);
  //   double y = r * sin(theta) * sin(phi);
  //   double z = r * cos(theta);
  //   cloud.insertNextPoint(hrleVectorType<double, D>(x, y, z));
  // }

  // random points in cube ----------------------------------------------------
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> dis(-5., 5.);
  // int numberOfPoints = 1000;
  // for(int i = 0; i < numberOfPoints; ++i) {
  //   double x = dis(gen);
  //   double y = dis(gen);
  //   double z = dis(gen);
  //   cloud.insertNextPoint(hrleVectorType<double, D>(x, y, z));
  // }

  // diamond -------------------------------------------------------------------
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, 0, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, 0, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, -1, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, 1, 0));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, 0, -1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(0, 0, 1));

  // cube ----------------------------------------------------------------------
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, -1, -2));
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, -1, 2));
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, 1, -2));
  // cloud.insertNextPoint(hrleVectorType<double, D>(-1, 1, 2));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, -1, -1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, -1, 1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, 1, -1));
  // cloud.insertNextPoint(hrleVectorType<double, D>(1, 1, 1));

  // cylinder ----------------------------------------------------------------
  unsigned numberOfBasePoints = 50;
  double radius = 5.;
  double height = 5.;
  for (unsigned i = 0; i < numberOfBasePoints; ++i) {
    double angle = 2 * PI * (double(i) / double(numberOfBasePoints));

    double x = radius * cos(angle);
    double y = radius * sin(angle);
    cloud->insertNextPoint(hrleVectorType<double, D>(x, y, -height));
    cloud->insertNextPoint(hrleVectorType<double, D>(x, y, height));
  }

  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsConvexHull<double, D>(mesh, cloud).apply();

  auto pointMesh = lsSmartPointer<lsMesh<>>::New();
  for (unsigned i = 0; i < cloud->points.size(); ++i) {
    pointMesh->nodes.push_back(std::array<double, 3>{
        cloud->points[i][0], cloud->points[i][1], cloud->points[i][2]});
    pointMesh->vertices.push_back(std::array<unsigned, 1>{i});
  }
  std::cout << "Output point cloud" << std::endl;
  lsVTKWriter<double>(pointMesh, lsFileFormatEnum::VTP, "points.vtp").apply();

  // std::cout << "Output surface mesh" << std::endl;
  // mesh.print();
  // lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP, "hull.vtp").apply();

  std::cout << "create level set" << std::endl;
  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New(0.18);
  lsMakeGeometry<double, D> geom;
  geom.setLevelSet(levelSet);
  geom.setGeometry(cloud);
  geom.apply();

  // now make into level set
  // lsDomain<double, D> levelSet(0.3);
  // lsFromSurfaceMesh<double, D>(levelSet, mesh).apply();

  auto LSMesh = lsSmartPointer<lsMesh<>>::New();
  std::cout << "Output level set grid" << std::endl;
  lsToMesh<double, D>(levelSet, LSMesh).apply();
  lsVTKWriter<double>(LSMesh, lsFileFormatEnum::VTP, "LS.vtp").apply();

  std::cout << "Output level set surface" << std::endl;
  lsToSurfaceMesh<double, D>(levelSet, LSMesh).apply();
  lsVTKWriter<double>(LSMesh, lsFileFormatEnum::VTP, "LSmesh.vtp").apply();

  // {
  //   lsMesh<> mesh;
  //   std::cout << "Extracting..." << std::endl;
  //   lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   mesh.print();
  //   lsVTKWriter<double>(mesh, "after2D.vtk").apply();
  // }

  return 0;
}
