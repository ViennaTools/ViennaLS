#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsReader.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriter.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

int main() {
  constexpr int D = 3;

  omp_set_num_threads(4);

  lsDomain<double, D> levelSet;

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, D>(levelSet, lsSphere<double, D>(centre, radius))
      .apply();

  lsPointData &data = levelSet.getPointData();
  typename lsPointData::ScalarDataType scalars;
  typename lsPointData::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet.getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename lsPointData::VectorDataType::value_type({double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "llaalalalalaalalalalalaalal");

  lsWriter<double, D>(levelSet, "test.lvst").apply();

  // read it in again
  lsDomain<double, D> newLevelSet;
  lsReader<double, D>(newLevelSet, "test.lvst").apply();

  lsMesh mesh;
  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "test.vtk").apply();

  return 0;
}
