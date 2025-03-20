#pragma once

#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// Abstract class defining the interface for
/// the velocity field used during advection using lsAdvect.
template <class T> class VelocityField {
public:
  VelocityField() = default;

  /// Should return a scalar value for the velocity at coordinate
  /// for a point of material with the given normalVector.
  virtual T getScalarVelocity(const Vec3D<T> & /*coordinate*/, int /*material*/,
                              const Vec3D<T> & /*normalVector*/,
                              unsigned long /*pointId*/) {
    return 0;
  }

  /// Like getScalarVelocity, but returns a velocity value for each
  /// cartesian direction.
  virtual Vec3D<T> getVectorVelocity(const Vec3D<T> & /*coordinate*/,
                                     int /*material*/,
                                     const Vec3D<T> & /*normalVector*/,
                                     unsigned long /*pointId*/) {
    return Vec3D<T>{0, 0, 0};
  }

  /// If lsLocalLaxFriedrichsAnalytical is used as the advection scheme,
  /// this is called to provide the analytical solution for the alpha
  /// values, needed for stable integration.
  virtual T getDissipationAlpha(int /*direction*/, int /*material*/,
                                const Vec3D<T> & /*centralDifferences*/) {
    return 0;
  }

  virtual ~VelocityField() = default;
};

} // namespace viennals
