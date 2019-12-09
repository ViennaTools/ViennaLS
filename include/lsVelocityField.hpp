#ifndef LS_VELOCITY_FIELD_HPP
#define LS_VELOCITY_FIELD_HPP

#include <array>

/// Abstract class defining the interface for
/// the velocity field used during advection using lsAdvect.
template <class T> class lsVelocityField {
public:
  lsVelocityField() {}

  virtual T getScalarVelocity(const std::array<T, 3> &coordinate, int material,
                              const std::array<T, 3> &normalVector) {
    return 0;
  }

  virtual std::array<T, 3>
  getVectorVelocity(const std::array<T, 3> &coordinate, int material,
                    const std::array<T, 3> &normalVector) {
    return {0, 0, 0};
  }

  virtual ~lsVelocityField() {}
};

#endif // LS_VELOCITY_FIELD_HPP
