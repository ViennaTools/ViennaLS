#ifndef LS_VELOCITY_FIELD_HPP
#define LS_VELOCITY_FIELD_HPP

#include <hrleVectorType.hpp>
#include <limits>

template <class T> class lsVelocityField {
public:
  virtual T getScalarVelocity(
      hrleVectorType<T, 3> coordinate, int material,
      hrleVectorType<T, 3> normalVector = hrleVectorType<T, 3>(0.)) = 0;

  virtual hrleVectorType<T, 3> getVectorVelocity(
      hrleVectorType<T, 3> coordinate, int material,
      hrleVectorType<T, 3> normalVector = hrleVectorType<T, 3>(0.)) = 0;
};

#endif // LS_VELOCITY_FIELD_HPP