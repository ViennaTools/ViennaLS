#ifndef LS_FAST_ADVECT_DISTRIBUTIONS_HPP
#define LS_FAST_ADVECT_DISTRIBUTIONS_HPP

#include <hrleVectorType.hpp>
#include <lsMessage.hpp>

/// Base class for distributions used by lsGeometricAdvect.
/// All functions are pure virtual and must be implemented
/// by any advection distribution.
template <class T, int D> class lsGeometricAdvectDistribution {
public:
  lsGeometricAdvectDistribution() {}

  /// Quick check whether a point relative to the distributions
  /// center is inside the distribution.
  virtual bool isInside(const std::array<hrleCoordType, 3> &initial,
                        const std::array<hrleCoordType, 3> &candidate,
                        double eps = 0.) const = 0;

  /// Returns the signed distance of a point relative to the distributions
  /// center. This is the signed manhatten distance to the nearest surface
  /// point.
  virtual T
  getSignedDistance(const std::array<hrleCoordType, 3> &initial,
                    const std::array<hrleCoordType, 3> &candidate) const = 0;

  /// Sets bounds to the bounding box of the distribution.
  virtual std::array<hrleCoordType, 6> getBounds() const = 0;

  virtual ~lsGeometricAdvectDistribution() {}
};

/// Concrete implementation of lsGeometricAdvectDistribution for a spherical
/// advection distribution.
template <class T, int D>
class lsSphereDistribution : public lsGeometricAdvectDistribution<T, D> {
public:
  const T radius = 0.;
  const T radius2;
  const T gridDelta;

  lsSphereDistribution(const T passedRadius, const T delta)
      : radius(passedRadius), radius2(radius * radius), gridDelta(delta) {}

  bool isInside(const std::array<hrleCoordType, 3> &initial,
                const std::array<hrleCoordType, 3> &candidate,
                double eps = 0.) const {
    hrleCoordType dot = 0.;
    for (unsigned i = 0; i < D; ++i) {
      double tmp = candidate[i] - initial[i];
      dot += tmp * tmp;
    }

    if (std::sqrt(dot) <= std::abs(radius) + eps)
      return true;
    else
      return false;
  }

  T getSignedDistance(const std::array<hrleCoordType, 3> &initial,
                      const std::array<hrleCoordType, 3> &candidate) const {
    T distance = std::numeric_limits<T>::max();
    std::array<hrleCoordType, D> v;
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
        if (D == 3)
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

  std::array<hrleCoordType, 6> getBounds() const {
    std::array<hrleCoordType, 6> bounds = {};
    for (unsigned i = 0; i < D; ++i) {
      bounds[2 * i] = -radius;
      bounds[2 * i + 1] = radius;
    }
    return bounds;
  }
};

/// Concrete implementation of lsGeometricAdvectDistribution
/// for a rectangular box distribution.
template <class T, int D>
class lsBoxDistribution : public lsGeometricAdvectDistribution<T, D> {
public:
  const hrleVectorType<T, 3> posExtent;
  const T gridDelta;

  lsBoxDistribution(const std::array<T, 3> &halfAxes, const T delta)
      : posExtent(halfAxes), gridDelta(delta) {
    for (unsigned i = 0; i < D; ++i) {
      if (posExtent[i] < gridDelta) {
        lsMessage::getInstance()
            .addWarning("One half-axis of lsBoxDistribution is smaller than "
                        "the grid Delta! This can lead to numerical errors "
                        "breaking the distribution!")
            .print();
      }
    }
  }

  bool isInside(const std::array<hrleCoordType, 3> &initial,
                const std::array<hrleCoordType, 3> &candidate,
                double eps = 0.) const {
    for (unsigned i = 0; i < D; ++i) {
      if (std::abs(candidate[i] - initial[i]) > (posExtent[i] + eps)) {
        return false;
      }
    }
    return true;
  }

  T getSignedDistance(const std::array<hrleCoordType, 3> &initial,
                      const std::array<hrleCoordType, 3> &candidate) const {
    T distance = std::numeric_limits<T>::lowest();
    for (unsigned i = 0; i < D; ++i) {
      T vector = std::abs(candidate[i] - initial[i]);
      distance = std::max(vector - posExtent[i], distance);
    }
    return distance;
  }

  std::array<hrleCoordType, 6> getBounds() const {
    std::array<hrleCoordType, 6> bounds = {};
    for (unsigned i = 0; i < D; ++i) {
      bounds[2 * i] = -posExtent[i];
      bounds[2 * i + 1] = posExtent[i];
    }
    return bounds;
  }
};

#endif // LS_FAST_ADVECT_DISTRIBUTIONS_HPP
