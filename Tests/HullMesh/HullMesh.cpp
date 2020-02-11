#include <iostream>

#include <lsConvexHull.hpp>
#include <lsGeometries.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 3;

  // omp_set_num_threads(4);

  lsPointCloud<double, D> cloud;

  cloud.insertNextPoint(hrleVectorType<double, D>(0, 0, 1));
  cloud.insertNextPoint(hrleVectorType<double, D>(1, 0, 1));
  cloud.insertNextPoint(hrleVectorType<double, D>(1, 1, 1));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, 1, 1));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, 0, 0.2));
  cloud.insertNextPoint(hrleVectorType<double, D>(1, 0, 0.2));
  cloud.insertNextPoint(hrleVectorType<double, D>(1, 1, 0.2));
  cloud.insertNextPoint(hrleVectorType<double, D>(0, 1, 0.2));

  lsMesh hull;
  lsConvexHull<double, D>(hull, cloud).apply();
  lsVTKWriter(hull, lsFileFormatEnum::VTP, "hull_" + std::to_string(1) + ".vtp")
      .apply();

  return 0;
}
