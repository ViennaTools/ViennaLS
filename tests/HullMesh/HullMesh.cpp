#include <iostream>

#include <lsConvexHull.hpp>
#include <lsGeometries.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  constexpr int D = 3;

  // omp_set_num_threads(4);

  auto cloud = ls::SmartPointer<ls::PointCloud<double, D>>::New();

  cloud->insertNextPoint({0, 0, 1});
  cloud->insertNextPoint({1, 0, 1});
  cloud->insertNextPoint({1, 1, 1});
  cloud->insertNextPoint({0, 1, 1});
  cloud->insertNextPoint({0, 0, 0.2});
  cloud->insertNextPoint({1, 0, 0.2});
  cloud->insertNextPoint({1, 1, 0.2});
  cloud->insertNextPoint({0, 1, 0.2});

  auto hull = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ConvexHull<double, D>(hull, cloud).apply();
  ls::VTKWriter<double>(hull, ls::FileFormatEnum::VTP,
                        "hull_" + std::to_string(1) + ".vtp")
      .apply();

  return 0;
}
