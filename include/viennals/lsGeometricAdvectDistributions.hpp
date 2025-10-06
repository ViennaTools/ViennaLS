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
                              unsigned long pointId) const = 0;

  /// Sets bounds to the bounding box of the distribution.
  virtual std::array<viennahrle::CoordType, 6> getBounds() const = 0;

  virtual bool useSurfacePointId() const { return false; }

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
                      unsigned long /*pointId*/) const override {
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

  bool useSurfacePointId() const override { return true; }
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
                      unsigned long /*pointId*/) const override {
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

  bool useSurfacePointId() const override { return true; }
};

template <class T, int D>
class CustomSphereDistribution : public GeometricAdvectDistribution<T, D> {

  const std::vector<T> radii_;
  const T gridDelta_;
  T maxRadius_ = 0;

public:
  CustomSphereDistribution(const std::vector<T> &radii, T delta)
      : radii_(radii), gridDelta_(delta) {
    for (unsigned i = 0; i < D; ++i) {
      maxRadius_ = std::max(maxRadius_, std::abs(radii_[i]));
    }
  }

  T getSignedDistance(const Vec3D<viennahrle::CoordType> &initial,
                      const Vec3D<viennahrle::CoordType> &candidate,
                      unsigned long pointId) const override {
    T distance = std::numeric_limits<T>::max();
    Vec3D<viennahrle::CoordType> v{};
    for (unsigned i = 0; i < D; ++i) {
      v[i] = candidate[i] - initial[i];
    }

    if (pointId >= radii_.size()) {
      Logger::getInstance()
          .addError("Point ID " + std::to_string(pointId) +
                    " is out of bounds for CustomSphereDistribution with " +
                    std::to_string(radii_.size()) + " radii!")
          .print();
      return distance;
    }
    const auto radius = radii_[pointId];
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
        T x = radius * radius - y * y - z * z;
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
      bounds[2 * i] = -maxRadius_;
      bounds[2 * i + 1] = maxRadius_;
    }
    return bounds;
  }

  bool useSurfacePointId() const override { return true; }
};

template <class T, int D>
class TrenchDistribution : public GeometricAdvectDistribution<T, D> {

  const T trenchWidth_;
  const T trenchDepth_;
  const T rate_;
  const T gridDelta_;

  const T bottomMed_;
  const T a_, b_, n_;

public:
  TrenchDistribution(const T trenchWidth, const T trenchDepth, const T rate,
                     const T gridDelta, const T bottomMed = 1.0,
                     const T a = 1.0, const T b = 1.0, const T n = 1.0)
      : trenchWidth_(trenchWidth), trenchDepth_(trenchDepth), rate_(rate),
        gridDelta_(gridDelta), bottomMed_(bottomMed), a_(a), b_(b), n_(n) {}

  T getSignedDistance(const Vec3D<viennahrle::CoordType> &initial,
                      const Vec3D<viennahrle::CoordType> &candidate,
                      unsigned long pointId) const override {
    T distance = std::numeric_limits<T>::max();
    Vec3D<viennahrle::CoordType> v{};
    for (unsigned i = 0; i < D; ++i) {
      v[i] = candidate[i] - initial[i];
    }

    T radius = 0;
    // if (std::abs(initial[D - 1]) < gridDelta_) {
    // radius = rate_;
    // } else
    if (std::abs(initial[D - 1] + trenchDepth_) < gridDelta_) {
      radius = bottomMed_;
    } else {
      radius =
          a_ * std::pow(1. - std::abs(initial[D - 1]) / trenchDepth_, n_) + b_;
      // radius = a_ * std::exp(-n_ * std::abs(initial[D - 1])) + b_;
    }

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
        T x = radius * radius - y * y - z * z;
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
      bounds[2 * i] = -rate_;
      bounds[2 * i + 1] = rate_;
    }
    return bounds;
  }

  bool useSurfacePointId() const override { return true; }
};
} // namespace viennals
