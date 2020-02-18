#ifndef LS_FAST_ADVECT_DISTRIBUTIONS_HPP
#define LS_FAST_ADVECT_DISTRIBUTIONS_HPP

/// Base class for distributions used by lsFastAdvect.
/// All functions are pure virtual and must be implemented
/// by any advection distribution.
template <class T, int D> class lsFastAdvectDistribution {
public:
  /// Quick check whether a point relative to the distributions
  /// center is inside the distribution.
  virtual bool isInside(hrleVectorType<T, D> &v, double eps = 0.) const = 0;

  /// Returns the signed distance of a point relative to the distributions
  /// center. This is the signed manhatten distance to the nearest surface point.
  virtual double getSignedDistance(hrleVectorType<T, D> &v) const = 0;

  /// Sets bounds to the bounding box of the distribution.
  virtual void getBounds(double *bounds) const = 0;
};

/// Concrete implementation of lsFastAdvectDistribution for a spherical
/// advection distribution.
template <class T, int D>
class lsSphereDistribution : public lsFastAdvectDistribution<T, D> {
public:
  const T radius = 0.;
  const T radius2;

  lsSphereDistribution(T passedRadius)
      : radius(passedRadius), radius2(radius * radius) {}

  bool isInside(hrleVectorType<T, D> &v, double eps = 0.) const {
    if (std::sqrt(DotProduct(v, v)) <= radius + eps)
      return true;
    else
      return false;
  }

  double getSignedDistance(hrleVectorType<T, D> &v) const {
    T distance = std::numeric_limits<T>::max();
    for (unsigned i = 0; i < D; ++i) {
      T y = (v[(i + 1) % D]);
      T z = 0;
      if (D == 3)
        z = (v[(i + 2) % D]);
      T x = radius2 - y * y - z * z;
      if (x < 0.)
        continue;
      T dirRadius =
          std::abs(v[i]) - std::sqrt(x);
      if (std::abs(dirRadius) < std::abs(distance))
        distance = dirRadius;
    }
    return distance;

    // return std::sqrt(DotProduct(v, v)) - radius;
  }

  void getBounds(double *bounds) const {
    for (unsigned i = 0; i < D; ++i) {
      bounds[2 * i] = -radius;
      bounds[2 * i + 1] = radius;
    }
  }
};

#endif // LS_FAST_ADVECT_DISTRIBUTIONS_HPP
