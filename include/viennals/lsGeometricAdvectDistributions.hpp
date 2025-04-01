#pragma once

#include <hrleTypes.hpp>

#include <vcLogger.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// Base class for distributions used by lsGeometricAdvect.
/// All functions are pure virtual and must be implemented
/// by any advection distribution.
template <class T, int D> class GeometricAdvectDistribution {
public:
  GeometricAdvectDistribution() = default;

  /// Quick check whether a point relative to the distributions
  /// center is inside the distribution. If there is no quick
  /// check due to the complexity of the distribution, always
  /// return true or do not overload this function.
  virtual bool isInside(const Vec3D<viennahrle::CoordType> &initial,
                        const Vec3D<viennahrle::CoordType> &candidate,
                        double eps) const {
    return true;
  }

  /// Returns the signed distance of a point relative to the distributions
  /// center. This is the signed manhatten distance to the nearest surface
  /// point.
  virtual T getSignedDistance(const Vec3D<viennahrle::CoordType> &initial,
                              const Vec3D<viennahrle::CoordType> &candidate,
                              unsigned long initialPointId) const = 0;

  /// Sets bounds to the bounding box of the distribution.
  virtual std::array<viennahrle::CoordType, 6> getBounds() const = 0;

  virtual ~GeometricAdvectDistribution() = default;
};

/// Concrete implementation of GeometricAdvectDistribution for a spherical
/// advection distribution.
template <class T, int D>
class SphereDistribution : public GeometricAdvectDistribution<T, D> {
public:
  const T radius = 0.;
  const T radius2;
  const T gridDelta;

  SphereDistribution(const T passedRadius, const T delta)
      : radius(passedRadius), radius2(radius * radius), gridDelta(delta) {}

  bool isInside(const Vec3D<viennahrle::CoordType> &initial,
                const Vec3D<viennahrle::CoordType> &candidate,
                double eps) const override {
    viennahrle::CoordType dot = 0.;
    for (unsigned i = 0; i < D; ++i) {
      double tmp = candidate[i] - initial[i];
      dot += tmp * tmp;
    }

    if (std::sqrt(dot) <= std::abs(radius) + eps)
      return true;
    else
      return false;
  }

  T getSignedDistance(const Vec3D<viennahrle::CoordType> &initial,
                      const Vec3D<viennahrle::CoordType> &candidate,
                      unsigned long /*initialPointId*/) const override {
    T distance = std::numeric_limits<T>::max();
    Vec3D<viennahrle::CoordType> v{};
    for (unsigned i = 0; i < D; ++i) {
      v[i] = candidate[i] - initial[i];
    }

    if (std::abs(radius) <= gridDelta) {
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

  std::array<viennahrle::CoordType, 6> getBounds() const override {
    std::array<viennahrle::CoordType, 6> bounds = {};
    for (unsigned i = 0; i < D; ++i) {
      bounds[2 * i] = -radius;
      bounds[2 * i + 1] = radius;
    }
    return bounds;
  }
};

/// Concrete implementation of GeometricAdvectDistribution
/// for a rectangular box distribution.
template <class T, int D>
class BoxDistribution : public GeometricAdvectDistribution<T, D> {
public:
  const VectorType<T, 3> posExtent;
  const T gridDelta;

  BoxDistribution(const std::array<T, 3> &halfAxes, const T delta)
      : posExtent(halfAxes), gridDelta(delta) {
    for (unsigned i = 0; i < D; ++i) {
      if (std::abs(posExtent[i]) < gridDelta) {
        Logger::getInstance()
            .addWarning("One half-axis of BoxDistribution is smaller than "
                        "the grid Delta! This can lead to numerical errors "
                        "breaking the distribution!")
            .print();
      }
    }
  }

  bool isInside(const Vec3D<viennahrle::CoordType> &initial,
                const Vec3D<viennahrle::CoordType> &candidate,
                double eps) const override {
    for (unsigned i = 0; i < D; ++i) {
      if (std::abs(candidate[i] - initial[i]) >
          (std::abs(posExtent[i]) + eps)) {
        return false;
      }
    }
    return true;
  }

  T getSignedDistance(const Vec3D<viennahrle::CoordType> &initial,
                      const Vec3D<viennahrle::CoordType> &candidate,
                      unsigned long /*initialPointId*/) const override {
    T distance = std::numeric_limits<T>::lowest();
    for (unsigned i = 0; i < D; ++i) {
      T vector = std::abs(candidate[i] - initial[i]);
      distance = std::max(vector - std::abs(posExtent[i]), distance);
    }
    return (posExtent[0] < 0) ? -distance : distance;
  }

  std::array<viennahrle::CoordType, 6> getBounds() const override {
    std::array<viennahrle::CoordType, 6> bounds = {};
    for (unsigned i = 0; i < D; ++i) {
      bounds[2 * i] = -posExtent[i];
      bounds[2 * i + 1] = posExtent[i];
    }
    return bounds;
  }
};

} // namespace viennals
