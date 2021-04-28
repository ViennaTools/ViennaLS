#include <iostream>

#include <lsConvexHull.hpp>
#include <lsGeometries.hpp>
#include <lsSmartPointer.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 3;

  // omp_set_num_threads(4);

  auto cloud = lsSmartPointer<lsPointCloud<double, D>>::New();

  cloud->insertNextPoint(hrleVectorType<double, D>(0, 0, 1));
  cloud->insertNextPoint(hrleVectorType<double, D>(1, 0, 1));
  cloud->insertNextPoint(hrleVectorType<double, D>(1, 1, 1));
  cloud->insertNextPoint(hrleVectorType<double, D>(0, 1, 1));
  cloud->insertNextPoint(hrleVectorType<double, D>(0, 0, 0.2));
  cloud->insertNextPoint(hrleVectorType<double, D>(1, 0, 0.2));
  cloud->insertNextPoint(hrleVectorType<double, D>(1, 1, 0.2));
  cloud->insertNextPoint(hrleVectorType<double, D>(0, 1, 0.2));

  auto hull = lsSmartPointer<lsMesh<>>::New();
  lsConvexHull<double, D>(hull, cloud).apply();
  lsVTKWriter<double>(hull, lsFileFormatEnum::VTP,
                      "hull_" + std::to_string(1) + ".vtp")
      .apply();

  return 0;
}
