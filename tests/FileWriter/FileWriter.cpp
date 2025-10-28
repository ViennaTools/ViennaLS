#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsReader.hpp>
#include <lsTestAsserts.hpp>
#include <lsWriter.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 10;
  double gridDelta = 1;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  const double radius = 7.3;
  const ls::VectorType<double, D> centre{5., 0.};
  // const double centre[D] = {0,0.1};
  // const double normal[D] = {0, 1};

  // ls::MakeGeometry<double, D>(levelSet, ls::SmartPointer<ls::Plane<double,
  // D>>::New(centre, normal)).apply();
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
      .apply();

  ls::PointData<double> &data = levelSet->getPointData();
  typename ls::PointData<double>::ScalarDataType scalars;
  typename ls::PointData<double>::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet->getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename ls::PointData<double>::VectorDataType::value_type(
            {double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "myVectors");

  // TODO: using the output of hrleDomain::print here to compare data in LS
  // there should definitely be a better way
  std::string domainString;
  {
    std::ostringstream oss;
    levelSet->getDomain().print(oss);
    domainString = oss.str();
  }

  ls::Writer<double, D>(levelSet, "test.lvst").apply();

  // read it in again
  auto newLevelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
  ls::Reader<double, D>(newLevelSet, "test.lvst").apply();

  {
    std::ostringstream oss;
    newLevelSet->getDomain().print(oss);
    VC_TEST_ASSERT(oss.str() == domainString)
  }

  return 0;
}
