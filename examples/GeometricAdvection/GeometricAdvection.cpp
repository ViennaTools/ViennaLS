#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  emulation, by creating a trench geometry. A uniform
  layer of a different material is then grown on top. It is
  the same example as Deposition but emulates the deposition
  rather than simulating a slow growth.
  \example GeometricAdvection.cpp
*/

namespace ls = viennals;

using NumericType = double;

template <class T, int D>
class MyDistribution : public ls::GeometricAdvectDistribution<T, D> {
private:
  const T gridDelta_ = 0.5;
  const T constantGrowth_ = 1;
  const T ARGrowth_ = 1;
  const T exponent_ = 1;
  const T height_;

public:
  MyDistribution(T gridDelta, T height, T constantGrowth, T ARGrowth,
                 T exponent)
      : gridDelta_(gridDelta), constantGrowth_(constantGrowth),
        ARGrowth_(ARGrowth), exponent_(exponent), height_(height) {}

  T getSignedDistance(const ls::Vec3D<viennahrle::CoordType> &initial,
                      const ls::Vec3D<viennahrle::CoordType> &candidate,
                      unsigned long pointId) const override {

    T distance = std::numeric_limits<T>::max();
    ls::Vec3D<viennahrle::CoordType> v{};
    for (unsigned i = 0; i < D; ++i) {
      v[i] = candidate[i] - initial[i];
    }

    const T radius =
        ARGrowth_ * std::pow(std::abs(initial[D - 1]) / height_, exponent_) +
        constantGrowth_;
    const T radius2 = radius * radius;
    if (std::abs(radius) <= gridDelta_) {
      distance =
          std::max(std::max(std::abs(v[0]), std::abs(v[1])), std::abs(v[2])) -
          std::abs(radius);
    } else {
      for (unsigned i = 0; i < D; ++i) {
        T y = (v[(i + 1) % D]);
        T z = 0;
        if constexpr (D == 3)
          z = (v[(i + 2) % D]);
        T x = radius2 - y * y - z * z;
        if (x < 0.)
          continue;
        T dirRadius = std::abs(v[i]) - std::sqrt(x);
        if (std::abs(dirRadius) < std::abs(distance))
          distance = dirRadius;
      }
    }
    // return distance;
    if (radius < 0) {
      return -distance;
    } else {
      return distance;
    }
  }

  /// Sets bounds to the bounding box of the distribution.
  std::array<viennahrle::CoordType, 6> getBounds() const override {
    return {{-(constantGrowth_ + ARGrowth_), (constantGrowth_ + ARGrowth_),
             -(constantGrowth_ + ARGrowth_), (constantGrowth_ + ARGrowth_),
             -(constantGrowth_ + ARGrowth_), (constantGrowth_ + ARGrowth_)}};
  }
};

int main() {

  constexpr int D = 2;
  omp_set_num_threads(1);

  NumericType extent = 30;
  NumericType gridDelta = 0.5;
  NumericType height = 50;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  {
    NumericType origin[D] = {0};
    origin[D - 1] = height;
    NumericType planeNormal[D] = {0.};
    planeNormal[D - 1] = 1.;
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    // make -x and +x greater than domain for numerical stability
    // NumericType ylimit = extent / 4.;
    NumericType minCorner[D] = {-extent / 3, 0.};
    NumericType maxCorner[D] = {extent / 3, height + 1.};
    auto box =
        ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<NumericType, D>(trench, box).apply();
    // Create trench geometry
    ls::BooleanOperation<NumericType, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "trench-0.vtp").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  auto newLayer = ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);

  std::cout << "Advecting" << std::endl;
  // Grow the layer uniformly by 4 as in deposition example
  auto dist = ls::SmartPointer<MyDistribution<viennahrle::CoordType, D>>::New(
      gridDelta, height, 2, 10, 3);
  ls::GeometricAdvect<NumericType, D>(newLayer, dist).apply();

  {
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "trench-final.vtp").apply();
  }

  return 0;
}
