#ifndef LS_VELOCITY_FIELD_HPP
#define LS_VELOCITY_FIELD_HPP

#include <hrleVectorType.hpp>
#include <limits>

/// Abstract class defining the interface for
/// the velocity field used during advection using lsAdvect.
template <class T> class lsVelocityField {
public:
  lsVelocityField() {}

  virtual T getScalarVelocity(
      hrleVectorType<T, 3> coordinate, int material,
      hrleVectorType<T, 3> normalVector = hrleVectorType<T, 3>(T(0))) = 0;

  virtual hrleVectorType<T, 3> getVectorVelocity(
      hrleVectorType<T, 3> coordinate, int material,
      hrleVectorType<T, 3> normalVector = hrleVectorType<T, 3>(T(0))) = 0;

  virtual ~lsVelocityField(){};
};

#endif // LS_VELOCITY_FIELD_HPP
