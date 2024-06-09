#ifndef LS_FINITE_DIFFERENCES_HPP
#define LS_FINITE_DIFFERENCES_HPP

#include <hrleVectorType.hpp>

#include <lsMessage.hpp>

namespace lsInternal {

// Numerical scheme definitions, central is default
enum class DifferentiationSchemeEnum : unsigned {
  FIRST_ORDER = 0,
  SECOND_ORDER = 1,
  WENO3 = 2,
  WENO5 = 3
};

template <class T, DifferentiationSchemeEnum scheme =
                       DifferentiationSchemeEnum::FIRST_ORDER>
class lsFiniteDifferences {
  template <class V> static V square(V x) { return x * x; }

public:
  lsFiniteDifferences(){};

  static unsigned getNumberOfValues(DifferentiationSchemeEnum s) {
    switch (s) {
    case DifferentiationSchemeEnum::FIRST_ORDER:
      return 3;
    case DifferentiationSchemeEnum::SECOND_ORDER:
    case DifferentiationSchemeEnum::WENO3:
      return 5;
    case DifferentiationSchemeEnum::WENO5:
      return 7;
    default:
      lsMessage::getInstance().addError("Invalid finite differences scheme!");
      return 0;
    }
  }

  /// Weighted essentially non-oscillatory differentiation scheme 3rd order
  /// x1 ... x5 stencil points from left to right
  /// plus == true => right-sided
  static T weno3(const T *x, T delta, bool plus, T eps = 1e-6) {
    T dx[4];
    for (unsigned i = 0; i < 4; ++i) {
      dx[i] = x[i + 1] - x[i];
    }

    T result = 0;
    if (plus) {
      T rp = (eps + lsFiniteDifferences::square(dx[3] - dx[2])) /
             (eps + lsFiniteDifferences::square(dx[2] - dx[1]));
      T wp = 1.0 / (1 + 2.0 * lsFiniteDifferences::square(rp));
      result = dx[1] + dx[2] - wp * (dx[3] - 2.0 * dx[2] + dx[1]);
    } else {
      T rp = (eps + lsFiniteDifferences::square(dx[1] - dx[0])) /
             (eps + lsFiniteDifferences::square(dx[2] - dx[1]));
      T wp = 1.0 / (1 + 2.0 * lsFiniteDifferences::square(rp));
      result = dx[1] + dx[2] - wp * (dx[0] - 2.0 * dx[1] + dx[2]);
    }

    return result / (2.0 * delta);
  }

  // Weighted essentially non-oscillatory differentiation scheme 5th order
  // x1 ... x7 stencil points from left to right
  // plus == true => right-sided
  static T weno5(const T *x, T dx, bool plus, T eps = 1e-6) {

    if (plus == false) {
      T v1 = (x[1] - x[0]) / dx; // i-3
      T v2 = (x[2] - x[1]) / dx; // i-2
      T v3 = (x[3] - x[2]) / dx; // i-1
      T v4 = (x[4] - x[3]) / dx; // i
      T v5 = (x[5] - x[4]) / dx; // i+1

      T p1 = v1 / 3.0 - 7 * v2 / 6.0 + 11 * v3 / 6.0;
      T p2 = -v2 / 6.0 + 5 * v3 / 6.0 + v4 / 3.0;
      T p3 = v3 / 3.0 + 5 * v4 / 6.0 - v5 / 6.0;

      T s1 = 13 / 12.0 * lsFiniteDifferences::square(v1 - 2 * v2 + v3) +
             1 / 4.0 * lsFiniteDifferences::square(v1 - 4 * v2 + 3 * v3);
      T s2 = 13 / 12.0 * lsFiniteDifferences::square(v2 - 2 * v3 + v4) +
             1 / 4.0 * lsFiniteDifferences::square(v2 - v4);
      T s3 = 13 / 12.0 * lsFiniteDifferences::square(v3 - 2 * v4 + v5) +
             1 / 4.0 * lsFiniteDifferences::square(3 * v3 - 4 * v4 + v5);

      T al1 = 0.1 / (eps + s1);
      T al2 = 0.6 / (eps + s2);
      T al3 = 0.3 / (eps + s3);

      T alsum = al1 + al2 + al3;

      T w1 = al1 / alsum;
      T w2 = al2 / alsum;
      T w3 = al3 / alsum;

      return w1 * p1 + w2 * p2 + w3 * p3;
    } else {
      T v1 = (x[6] - x[5]) / dx;
      T v2 = (x[5] - x[4]) / dx;
      T v3 = (x[4] - x[3]) / dx;
      T v4 = (x[3] - x[2]) / dx;
      T v5 = (x[2] - x[1]) / dx;

      T p1 = v1 / 3.0 - 7 * v2 / 6.0 + 11 * v3 / 6.0;
      T p2 = -v2 / 6.0 + 5 * v3 / 6.0 + v4 / 3.0;
      T p3 = v3 / 3.0 + 5 * v4 / 6.0 - v5 / 6.0;

      T s1 = 13 / 12.0 * lsFiniteDifferences::square(v1 - 2 * v2 + v3) +
             1 / 4.0 * lsFiniteDifferences::square(v1 - 4 * v2 + 3 * v3);
      T s2 = 13 / 12.0 * lsFiniteDifferences::square(v2 - 2 * v3 + v4) +
             1 / 4.0 * lsFiniteDifferences::square(v2 - v4);
      T s3 = 13 / 12.0 * lsFiniteDifferences::square(v3 - 2 * v4 + v5) +
             1 / 4.0 * lsFiniteDifferences::square(3 * v3 - 4 * v4 + v5);

      T al1 = 0.1 / (eps + s1);
      T al2 = 0.6 / (eps + s2);
      T al3 = 0.3 / (eps + s3);

      T alsum = al1 + al2 + al3;

      T w1 = al1 / alsum;
      T w2 = al2 / alsum;
      T w3 = al3 / alsum;

      return w1 * p1 + w2 * p2 + w3 * p3;
    }
  }

  // Finite difference in the negative direction using the scheme specified
  // by scheme. The passed vector contains the required neighbouring values,
  // with the center point being in the centre of the vector.
  // First order: x_-1, x, x_+1
  // Second order: x_-2, x_-1, x, x_+1, x_+2
  // WENO3: x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  static T differenceNegative(const T *values, const double &delta) {
    if constexpr (scheme == DifferentiationSchemeEnum::FIRST_ORDER) {
      return (values[1] - values[0]) / delta;
    } else if (scheme == DifferentiationSchemeEnum::SECOND_ORDER) {
      // TODO: implement second order integration here
    } else if (scheme == DifferentiationSchemeEnum::WENO3) {
      return weno3(values, delta, false);
    } else if (scheme == DifferentiationSchemeEnum::WENO5) {
      return weno5(values, delta, false);
    } else {
      lsMessage::getInstance().addError("Invalid finite differences scheme!");
    }
  }

  // Finite difference in the negative direction using the scheme specified
  // by scheme. The passed vector contains the required neighbouring values,
  // with the center point being in the centre of the vector.
  // First order:       x_-1, x, x_+1
  // WENO3:       x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  static T differencePositive(const T *values, const double &delta) {
    if constexpr (scheme == DifferentiationSchemeEnum::FIRST_ORDER) {
      return (values[2] - values[1]) / delta;
    } else if (scheme == DifferentiationSchemeEnum::SECOND_ORDER) {
      // TODO: implement second order integration here
    } else if (scheme == DifferentiationSchemeEnum::WENO3) {
      return weno3(values, delta, true);
    } else if (scheme == DifferentiationSchemeEnum::WENO5) {
      return weno5(values, delta, true);
    } else {
      lsMessage::getInstance().addError("Invalid finite differences scheme!");
    }
  }

  // Calculates the gradient around the central point in the passed vector
  // Depending on the scheme, the passed vector takes a different size:
  // First order:       x_-1, x, x_+1
  // WENO3:       x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  static T calculateGradient(const T *values, const double &delta) {
    return (differencePositive(values, delta) +
            differenceNegative(values, delta)) *
           0.5;
  }

  static T calculateGradientDiff(const T *values, const double &delta) {
    return (differencePositive(values, delta) -
            differenceNegative(values, delta)) *
           0.5;
  }
};

} // namespace lsInternal

#endif // LS_FINITE_DIFFERENCES_HPP
