#pragma once

#include <vcLogger.hpp>

namespace lsInternal {

using namespace viennacore;

// Numerical scheme definitions, central is default
enum class DifferentiationSchemeEnum : unsigned {
  FIRST_ORDER = 0,
  SECOND_ORDER = 1,
  WENO3 = 2,
  WENO5 = 3,
  WENO5_Z = 4
};

template <class T, DifferentiationSchemeEnum scheme =
                       DifferentiationSchemeEnum::FIRST_ORDER>
class FiniteDifferences {
  template <class V> static V square(V x) { return x * x; }

public:
  FiniteDifferences() = default;

  static unsigned getNumberOfValues(DifferentiationSchemeEnum s) {
    switch (s) {
    case DifferentiationSchemeEnum::FIRST_ORDER:
      return 3;
    case DifferentiationSchemeEnum::SECOND_ORDER:
    case DifferentiationSchemeEnum::WENO3:
      return 5;
    case DifferentiationSchemeEnum::WENO5:
    case DifferentiationSchemeEnum::WENO5_Z:
      return 7;
    default:
      Logger::getInstance().addError("Invalid finite differences scheme!");
      return 0;
    }
  }

  static constexpr unsigned getNumberOfValues() {
    switch (scheme) {
    case DifferentiationSchemeEnum::FIRST_ORDER:
      return 3;
    case DifferentiationSchemeEnum::SECOND_ORDER:
    case DifferentiationSchemeEnum::WENO3:
      return 5;
    case DifferentiationSchemeEnum::WENO5:
    case DifferentiationSchemeEnum::WENO5_Z:
      return 7;
    default:
      Logger::getInstance().addError("Invalid finite differences scheme!");
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
      T num = eps + FiniteDifferences::square(dx[3] - dx[2]);
      T den = eps + FiniteDifferences::square(dx[2] - dx[1]);
      T wp = (den * den) / (den * den + 2.0 * num * num);
      result = dx[1] + dx[2] - wp * (dx[3] - 2.0 * dx[2] + dx[1]);
    } else {
      T num = eps + FiniteDifferences::square(dx[1] - dx[0]);
      T den = eps + FiniteDifferences::square(dx[2] - dx[1]);
      T wp = (den * den) / (den * den + 2.0 * num * num);
      result = dx[1] + dx[2] - wp * (dx[0] - 2.0 * dx[1] + dx[2]);
    }

    return result / (2.0 * delta);
  }

  // Weighted essentially non-oscillatory differentiation scheme 5th order
  // x1 ... x7 stencil points from left to right
  // plus == true => right-sided
  static T weno5(const T *x, T dx, bool plus, T eps = 1e-6) {
    const T inv_dx = T(1) / dx;
    const T c1_6 = T(1) / T(6);
    const T c13_12 = T(13) / T(12);
    const T c1_4 = T(1) / T(4);

    if (plus == false) {
      T v1 = (x[1] - x[0]) * inv_dx; // i-3
      T v2 = (x[2] - x[1]) * inv_dx; // i-2
      T v3 = (x[3] - x[2]) * inv_dx; // i-1
      T v4 = (x[4] - x[3]) * inv_dx; // i
      T v5 = (x[5] - x[4]) * inv_dx; // i+1

      T p1 = (T(2) * v1 - T(7) * v2 + T(11) * v3) * c1_6;
      T p2 = (-v2 + T(5) * v3 + T(2) * v4) * c1_6;
      T p3 = (T(2) * v3 + T(5) * v4 - v5) * c1_6;

      T s1 = c13_12 * FiniteDifferences::square(v1 - T(2) * v2 + v3) +
             c1_4 * FiniteDifferences::square(v1 - T(4) * v2 + T(3) * v3);
      T s2 = c13_12 * FiniteDifferences::square(v2 - T(2) * v3 + v4) +
             c1_4 * FiniteDifferences::square(v2 - v4);
      T s3 = c13_12 * FiniteDifferences::square(v3 - T(2) * v4 + v5) +
             c1_4 * FiniteDifferences::square(T(3) * v3 - T(4) * v4 + v5);

      T al1 = T(0.1) / (eps + s1);
      T al2 = T(0.6) / (eps + s2);
      T al3 = T(0.3) / (eps + s3);

      return (al1 * p1 + al2 * p2 + al3 * p3) / (al1 + al2 + al3);
    } else {
      T v1 = (x[6] - x[5]) * inv_dx;
      T v2 = (x[5] - x[4]) * inv_dx;
      T v3 = (x[4] - x[3]) * inv_dx;
      T v4 = (x[3] - x[2]) * inv_dx;
      T v5 = (x[2] - x[1]) * inv_dx;

      T p1 = (T(2) * v1 - T(7) * v2 + T(11) * v3) * c1_6;
      T p2 = (-v2 + T(5) * v3 + T(2) * v4) * c1_6;
      T p3 = (T(2) * v3 + T(5) * v4 - v5) * c1_6;

      T s1 = c13_12 * FiniteDifferences::square(v1 - T(2) * v2 + v3) +
             c1_4 * FiniteDifferences::square(v1 - T(4) * v2 + T(3) * v3);
      T s2 = c13_12 * FiniteDifferences::square(v2 - T(2) * v3 + v4) +
             c1_4 * FiniteDifferences::square(v2 - v4);
      T s3 = c13_12 * FiniteDifferences::square(v3 - T(2) * v4 + v5) +
             c1_4 * FiniteDifferences::square(T(3) * v3 - T(4) * v4 + v5);

      T al1 = T(0.1) / (eps + s1);
      T al2 = T(0.6) / (eps + s2);
      T al3 = T(0.3) / (eps + s3);

      return (al1 * p1 + al2 * p2 + al3 * p3) / (al1 + al2 + al3);
    }
  }

  // Weighted essentially non-oscillatory differentiation scheme 5th order (Z)
  // x1 ... x7 stencil points from left to right
  // plus == true => right-sided
  static T weno5z(const T *x, T dx, bool plus, T eps = 1e-6) {
    const T inv_dx = T(1) / dx;
    const T c1_6 = T(1) / T(6);
    const T c13_12 = T(13) / T(12);
    const T c1_4 = T(1) / T(4);

    if (plus == false) {
      T v1 = (x[1] - x[0]) * inv_dx; // i-3
      T v2 = (x[2] - x[1]) * inv_dx; // i-2
      T v3 = (x[3] - x[2]) * inv_dx; // i-1
      T v4 = (x[4] - x[3]) * inv_dx; // i
      T v5 = (x[5] - x[4]) * inv_dx; // i+1

      T p1 = (T(2) * v1 - T(7) * v2 + T(11) * v3) * c1_6;
      T p2 = (-v2 + T(5) * v3 + T(2) * v4) * c1_6;
      T p3 = (T(2) * v3 + T(5) * v4 - v5) * c1_6;

      T s1 = c13_12 * FiniteDifferences::square(v1 - T(2) * v2 + v3) +
             c1_4 * FiniteDifferences::square(v1 - T(4) * v2 + T(3) * v3);
      T s2 = c13_12 * FiniteDifferences::square(v2 - T(2) * v3 + v4) +
             c1_4 * FiniteDifferences::square(v2 - v4);
      T s3 = c13_12 * FiniteDifferences::square(v3 - T(2) * v4 + v5) +
             c1_4 * FiniteDifferences::square(T(3) * v3 - T(4) * v4 + v5);

      T tau5 = std::abs(s1 - s3);
      T tau5sq = FiniteDifferences::square(tau5);

      T al1 = T(0.1) * (T(1) + tau5sq / FiniteDifferences::square(eps + s1));
      T al2 = T(0.6) * (T(1) + tau5sq / FiniteDifferences::square(eps + s2));
      T al3 = T(0.3) * (T(1) + tau5sq / FiniteDifferences::square(eps + s3));

      return (al1 * p1 + al2 * p2 + al3 * p3) / (al1 + al2 + al3);
    } else {
      T v1 = (x[6] - x[5]) * inv_dx;
      T v2 = (x[5] - x[4]) * inv_dx;
      T v3 = (x[4] - x[3]) * inv_dx;
      T v4 = (x[3] - x[2]) * inv_dx;
      T v5 = (x[2] - x[1]) * inv_dx;

      T p1 = (T(2) * v1 - T(7) * v2 + T(11) * v3) * c1_6;
      T p2 = (-v2 + T(5) * v3 + T(2) * v4) * c1_6;
      T p3 = (T(2) * v3 + T(5) * v4 - v5) * c1_6;

      T s1 = c13_12 * FiniteDifferences::square(v1 - T(2) * v2 + v3) +
             c1_4 * FiniteDifferences::square(v1 - T(4) * v2 + T(3) * v3);
      T s2 = c13_12 * FiniteDifferences::square(v2 - T(2) * v3 + v4) +
             c1_4 * FiniteDifferences::square(v2 - v4);
      T s3 = c13_12 * FiniteDifferences::square(v3 - T(2) * v4 + v5) +
             c1_4 * FiniteDifferences::square(T(3) * v3 - T(4) * v4 + v5);

      T tau5 = std::abs(s1 - s3);
      T tau5sq = FiniteDifferences::square(tau5);

      T al1 = T(0.1) * (T(1) + tau5sq / FiniteDifferences::square(eps + s1));
      T al2 = T(0.6) * (T(1) + tau5sq / FiniteDifferences::square(eps + s2));
      T al3 = T(0.3) * (T(1) + tau5sq / FiniteDifferences::square(eps + s3));

      return (al1 * p1 + al2 * p2 + al3 * p3) / (al1 + al2 + al3);
    }
  }

  // Finite difference in the negative direction using the scheme specified
  // by scheme. The passed vector contains the required neighbouring values,
  // with the center point being in the centre of the vector.
  // First order: x_-1, x, x_+1
  // Second order: x_-2, x_-1, x, x_+1, x_+2
  // WENO3: x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  // WENO5_Z: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  static T differenceNegative(const T *values, const double &delta) {
    if constexpr (scheme == DifferentiationSchemeEnum::FIRST_ORDER) {
      return (values[1] - values[0]) / delta;
    } else if (scheme == DifferentiationSchemeEnum::SECOND_ORDER) {
      // TODO: implement second order differentiation here
      Logger::getInstance().addError("Second order scheme not implemented!");
      return 0;
    } else if (scheme == DifferentiationSchemeEnum::WENO3) {
      return weno3(values, delta, false);
    } else if (scheme == DifferentiationSchemeEnum::WENO5) {
      return weno5(values, delta, false);
    } else if (scheme == DifferentiationSchemeEnum::WENO5_Z) {
      return weno5z(values, delta, false);
    } else {
      Logger::getInstance().addError("Invalid finite differences scheme!");
      return 0;
    }
  }

  // Finite difference in the negative direction using the scheme specified
  // by scheme. The passed vector contains the required neighbouring values,
  // with the center point being in the centre of the vector.
  // First order:       x_-1, x, x_+1
  // WENO3:       x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  // WENO5_Z: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  static T differencePositive(const T *values, const double &delta) {
    if constexpr (scheme == DifferentiationSchemeEnum::FIRST_ORDER) {
      return (values[2] - values[1]) / delta;
    } else if (scheme == DifferentiationSchemeEnum::SECOND_ORDER) {
      // TODO: implement second order differentiation here
      Logger::getInstance().addError("Second order scheme not implemented!");
      return 0;
    } else if (scheme == DifferentiationSchemeEnum::WENO3) {
      return weno3(values, delta, true);
    } else if (scheme == DifferentiationSchemeEnum::WENO5) {
      return weno5(values, delta, true);
    } else if (scheme == DifferentiationSchemeEnum::WENO5_Z) {
      return weno5z(values, delta, true);
    } else {
      Logger::getInstance().addError("Invalid finite differences scheme!");
      return 0;
    }
  }

  // Calculates the gradient around the central point in the passed vector
  // Depending on the scheme, the passed vector takes a different size:
  // First order:       x_-1, x, x_+1
  // WENO3:       x_-2, x_-1, x, x_+1, x_+2
  // WENO5: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
  // WENO5_Z: x_-3, x_-2, x_-1, x, x_+1, x_+2, x_+3
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
