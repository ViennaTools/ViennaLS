#pragma once

#include <vcLogger.hpp>

namespace lsInternal {

using namespace viennacore;

// Numerical scheme definitions, central is default
enum class DifferentiationSchemeEnum : unsigned {
  FIRST_ORDER = 0,
  SECOND_ORDER = 1,
  WENO3 = 2,
  WENO5 = 3
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
    T d[3];
    if (plus) {
      d[0] = x[2] - x[1];
      d[1] = x[3] - x[2];
      d[2] = x[4] - x[3];

      T N = eps + FiniteDifferences::square(d[2] - d[1]);
      T D = eps + FiniteDifferences::square(d[1] - d[0]);
      T D2 = D * D;
      T wp = D2 / (D2 + 2 * N * N);

      return (d[0] + d[1] - wp * (d[2] - 2 * d[1] + d[0])) / (2 * delta);
    } else {
      d[0] = x[1] - x[0];
      d[1] = x[2] - x[1];
      d[2] = x[3] - x[2];

      T N = eps + FiniteDifferences::square(d[1] - d[0]);
      T D = eps + FiniteDifferences::square(d[2] - d[1]);
      T D2 = D * D;
      T wp = D2 / (D2 + 2 * N * N);

      return (d[1] + d[2] - wp * (d[0] - 2 * d[1] + d[2])) / (2 * delta);
    }
  }

  // Weighted essentially non-oscillatory differentiation scheme 5th order
  // x1 ... x7 stencil points from left to right
  // plus == true => right-sided
  static T weno5(const T *x, T dx, bool plus, T eps = 1e-6) {
    // Optimized implementation avoiding multiple divisions by dx
    T d[5];
    if (!plus) {
      for (int i = 0; i < 5; ++i)
        d[i] = x[i + 1] - x[i];
    } else {
      for (int i = 0; i < 5; ++i)
        d[i] = x[6 - i] - x[5 - i];
    }

    // Smoothness indicators (scaled by dx^2)
    T s1 = T(13.0 / 12.0) * FiniteDifferences::square(d[0] - 2 * d[1] + d[2]) +
           T(1.0 / 4.0) * FiniteDifferences::square(d[0] - 4 * d[1] + 3 * d[2]);
    T s2 = T(13.0 / 12.0) * FiniteDifferences::square(d[1] - 2 * d[2] + d[3]) +
           T(1.0 / 4.0) * FiniteDifferences::square(d[1] - d[3]);
    T s3 = T(13.0 / 12.0) * FiniteDifferences::square(d[2] - 2 * d[3] + d[4]) +
           T(1.0 / 4.0) * FiniteDifferences::square(3 * d[2] - 4 * d[3] + d[4]);

    // Scale epsilon by dx^2 to match scaled smoothness indicators
    T eps_scaled = eps * dx * dx;

    T al1 = T(0.1) / (eps_scaled + s1);
    T al2 = T(0.6) / (eps_scaled + s2);
    T al3 = T(0.3) / (eps_scaled + s3);

    T alsum = al1 + al2 + al3;

    // Polynomials (scaled by dx)
    T p1 = T(1.0 / 6.0) * (2 * d[0] - 7 * d[1] + 11 * d[2]);
    T p2 = T(1.0 / 6.0) * (-d[1] + 5 * d[2] + 2 * d[3]);
    T p3 = T(1.0 / 6.0) * (2 * d[2] + 5 * d[3] - d[4]);

    // Final result divides by dx once
    return (al1 * p1 + al2 * p2 + al3 * p3) / (alsum * dx);
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
      // TODO: implement second order differentiation here
      Logger::getInstance().addError("Second order scheme not implemented!");
      return 0;
    } else if (scheme == DifferentiationSchemeEnum::WENO3) {
      return weno3(values, delta, false);
    } else if (scheme == DifferentiationSchemeEnum::WENO5) {
      return weno5(values, delta, false);
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
