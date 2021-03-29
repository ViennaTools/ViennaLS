#ifndef LS_VELOCITY_FIELD_HPP
#define LS_VELOCITY_FIELD_HPP

#include <array>

/// Abstract class defining the interface for
/// the velocity field used during advection using lsAdvect.
template <class T> class lsVelocityField {
public:
  lsVelocityField() {}

  /// Should return a scalar value for the velocity at coordinate
  /// for a point of material with the given normalVector.
  virtual T getScalarVelocity(const std::array<T, 3> & /*coordinate*/,
                              int /*material*/,
                              const std::array<T, 3> & /*normalVector*/,
                              unsigned long /*pointId*/) {
    return 0;
  }

  /// Like getScalarVelocity, but returns a velocity value for each
  /// cartesian direction.
  virtual std::array<T, 3>
  getVectorVelocity(const std::array<T, 3> & /*coordinate*/, int /*material*/,
                    const std::array<T, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return {0, 0, 0};
  }

  /// If lsLocalLaxFriedrichsAnalytical is used as the advection scheme,
  /// this is called to provide the analytical solution for the alpha
  /// values, needed for stable integration.
  virtual T
  getDissipationAlpha(int /*direction*/, int /*material*/,
                      const std::array<T, 3> & /*centralDifferences*/) {
    return 0;
  }

  virtual ~lsVelocityField() {}
};

#endif // LS_VELOCITY_FIELD_HPP
