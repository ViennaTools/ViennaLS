#pragma once

#include <array>

#include <vcVectorUtil.hpp>

namespace viennals {

using namespace viennacore;

/// Abstract class defining the interface for
/// the velocity field used during advection using lsAdvect.
template <class T> class VelocityField {
public:
  VelocityField() {}

  /// Should return a scalar value for the velocity at coordinate
  /// for a point of material with the given normalVector.
  virtual T getScalarVelocity(const Triple<T> & /*coordinate*/,
                              int /*material*/,
                              const Triple<T> & /*normalVector*/,
                              unsigned long /*pointId*/) {
    return 0;
  }

  /// Like getScalarVelocity, but returns a velocity value for each
  /// cartesian direction.
  virtual Triple<T> getVectorVelocity(const Triple<T> & /*coordinate*/,
                                      int /*material*/,
                                      const Triple<T> & /*normalVector*/,
                                      unsigned long /*pointId*/) {
    return {0, 0, 0};
  }

  /// If lsLocalLaxFriedrichsAnalytical is used as the advection scheme,
  /// this is called to provide the analytical solution for the alpha
  /// values, needed for stable integration.
  virtual T getDissipationAlpha(int /*direction*/, int /*material*/,
                                const Triple<T> & /*centralDifferences*/) {
    return 0;
  }

  virtual ~VelocityField() {}
};

} // namespace viennals
