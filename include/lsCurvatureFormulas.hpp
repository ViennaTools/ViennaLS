#ifndef LS_CURVATURE_FORMULAS_HPP
#define LS_CURVATURE_FORMULAS_HPP

#include <cmath>

#include <lsDomain.hpp>

namespace lsInternal {

// Formulas for space curves and higher dimension can be found here:
// https://doi.org/10.1016/j.cagd.2005.06.005

/// Returns the squared sum square for values in the range [start, end[
template <class It> double squareSumSquare(It begin, It end) {
  double result = 0.;
  while (begin != end) {
    result += *begin * *begin;
    ++begin;
  }
  return result * result;
}

/// Returns the root sum square to the 3rd power for values in the range [start,
/// end[
template <class It> double rootSumSquarePow3(It begin, It end) {
  double result = 0.;
  while (begin != end) {
    result += *begin * *begin;
    ++begin;
  }
  return std::pow(result, 1.5);
}

/// Mean curvature formula for implicit surfaces in 2D, the passed array should
/// contain the function values in the following order: (F_x, F_y, F_z, F_xx,
/// F_yy, F_zz, F_xy, F_yz, F_zx)
template <class T, std::size_t N>
double meanCurvature2D(std::array<T, N> funcValues) {
  double norm = rootSumSquarePow3(funcValues.begin(), funcValues.begin() + 2);
  auto &d = funcValues;

  return (d[3] * d[1] * d[1] - 2. * d[1] * d[0] * d[6] + d[4] * d[0] * d[0]) /
         norm;
}

/// Mean curvature formula for implicit surfaces in 3D, the passed array should
/// contain the function values in the following order: (F_x, F_y, F_z, F_xx,
/// F_yy, F_zz, F_xy, F_yz, F_zx)
template <class T, std::size_t N>
double meanCurvature3D(std::array<T, N> funcValues) {
  double norm = rootSumSquarePow3(funcValues.begin(), funcValues.begin() + 3);
  auto &d = funcValues;

  return (d[0] * d[0] * (d[4] + d[5]) + d[1] * d[1] * (d[3] + d[5]) +
          d[2] * d[2] * (d[3] + d[4]) +
          -2. *
              (d[0] * d[1] * d[6] + d[0] * d[2] * d[8] + d[1] * d[2] * d[7])) /
         (2. * norm);
}

/// Gaussian curvature formula for implicit surfaces in 3D, the passed array
/// should contain the function values in the following order: (F_x, F_y, F_z,
/// F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
template <class T, std::size_t N>
double gaussianCurvature3D(std::array<T, N> funcValues) {
  double norm = squareSumSquare(funcValues.begin(), funcValues.begin() + 3);
  auto &d = funcValues;

  return -(d[0] * d[0] * (d[7] * d[7] - d[4] * d[5]) +
           d[1] * d[1] * (d[8] * d[8] - d[3] * d[5]) +
           d[2] * d[2] * (d[6] * d[6] - d[3] * d[4]) +
           2. * (d[0] * d[1] * (d[5] * d[6] - d[8] * d[7]) +
                 d[0] * d[2] * (d[4] * d[8] - d[6] * d[7]) +
                 d[1] * d[2] * (d[3] * d[7] - d[6] * d[8]))) /
         norm;
}

/// Fills an std::array with differential values calculated from
/// neighbour values. This stencil only uses direct neighbours
/// for fast calculation of the differentials.
template <class It, class T = typename It::DomainType::hrleValueType>
std::array<T, 9> smallStencilFromIterator(It &it, const double gridDelta) {
  constexpr int D = It::DomainType::dimension;
  std::array<T, 9> d;
  for (int i = 0; i < D; i++) {
    hrleVectorType<hrleIndexType, D> posUnit(0);
    hrleVectorType<hrleIndexType, D> negUnit(0);
    posUnit[i] = 1;
    negUnit[i] = -1;
    T phi_0 = it.getCenter().getValue();
    T phi_px = it.getNeighbor(posUnit).getValue();
    T phi_nx = it.getNeighbor(negUnit).getValue();
    int second_pos = (i + 1) % D;
    posUnit[second_pos] = 1;
    negUnit[second_pos] = 1;
    T phi_pp = it.getNeighbor(posUnit).getValue();
    T phi_np = it.getNeighbor(negUnit).getValue();
    posUnit[second_pos] = -1;
    negUnit[second_pos] = -1;
    T phi_pn = it.getNeighbor(posUnit).getValue();
    T phi_nn = it.getNeighbor(negUnit).getValue();

    // Calculate derivatives for Hessian
    d[i] = (phi_px - phi_nx) / 2.;
    d[i + 3] = (phi_px - 2. * phi_0 + phi_nx) / gridDelta;
    d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) / (4 * gridDelta);
  }
  return std::move(d);
}

/// Fills an std::array with differential values calculated from
/// neighbour values. This stencil also uses diagonal neighbours
/// to achieve a higher accuracy.
template <class It, class T = typename It::DomainType::hrleValueType>
std::array<T, 9> bigStencilFromIterator(It &it, const double gridDelta) {
  constexpr int D = It::DomainType::dimension;
  std::array<T, 9> d;
  const double gridDelta2 = gridDelta * gridDelta;
  for (int i = 0; i < D; i++) {
    hrleVectorType<hrleIndexType, D> posUnit(0);
    hrleVectorType<hrleIndexType, D> negUnit(0);
    int first_axis = i;
    int second_axis = (i + 1) % D;
    posUnit[first_axis] = 1;
    negUnit[first_axis] = -1;
    T phi_0 = it.getCenter().getValue();
    T phi_px = it.getNeighbor(posUnit).getValue();
    T phi_nx = it.getNeighbor(negUnit).getValue();
    posUnit[second_axis] = 1;
    negUnit[second_axis] = 1;
    T phi_pp = it.getNeighbor(posUnit).getValue();
    T phi_np = it.getNeighbor(negUnit).getValue();
    posUnit[second_axis] = -1;
    negUnit[second_axis] = -1;
    T phi_pn = it.getNeighbor(posUnit).getValue();
    T phi_nn = it.getNeighbor(negUnit).getValue();
    posUnit[first_axis] = 0;
    negUnit[first_axis] = 0;
    posUnit[second_axis] = 1;
    negUnit[second_axis] = -1;
    T phi_py = it.getNeighbor(posUnit).getValue();
    T phi_ny = it.getNeighbor(negUnit).getValue();

    d[i] = (phi_pp - phi_np + phi_pn - phi_nn) / (4 * gridDelta);
    d[i + 3] = (phi_pp - 2. * phi_py + phi_np + phi_px - 2. * phi_0 + phi_nx +
                phi_pn - 2. * phi_ny + phi_nn) /
               (3 * gridDelta2);
    d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) / (4 * gridDelta2);
  }
  return std::move(d);
}

/// Calculates the Mean Curvature and/or the Gaussian Curvature of the level set
/// function in 3D Calculates the Curvature of the level set function in 2D.
/// Requires an iterator that is big enough to calculate second order
/// derivatives(e.g. hrleBoxIterator or hrleCartesianPlaneIterator)
template <class T, int D> class curvatureGeneralFormula {
private:
  T gridDelta;

public:
  curvatureGeneralFormula(T mGD) : gridDelta(mGD) {}

  template <class IteratorType>
  T meanCurvature(IteratorType &neighborIterator) {
    std::array<T, 9> d = smallStencilFromIterator(neighborIterator, gridDelta);
    if constexpr (D == 2) {
      return meanCurvature2D(d);
    } else {
      return meanCurvature3D(d);
    }
  }

  template <class IteratorType>
  T gaussianCurvature(IteratorType &neighborIterator) {
    if constexpr (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }
    std::array<T, 9> d = smallStencilFromIterator(neighborIterator, gridDelta);
    return gaussianCurvature3D(d);
  }

  template <class IteratorType>
  std::array<T, 2> meanGaussianCurvature(IteratorType &neighborIterator) {
    if constexpr (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }
    std::array<T, 9> d = smallStencilFromIterator(neighborIterator, gridDelta);
    T mean = meanCurvature3D(d);
    T gauss = gaussianCurvature3D(d);

    return std::array<T, 2>{mean, gauss};
  }
};

// Calculates the Curvature using the General Formula for implict surfaces and
// uses different first order approximations for D_x and D_xx derivatives
template <class T, int D> class curvatureGeneralFormulaBigStencil {

private:
  T gridDelta;

public:
  curvatureGeneralFormulaBigStencil(T mGD) : gridDelta(mGD) {}

  template <class IteratorType>
  T meanCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d = bigStencilFromIterator(neighborIterator, gridDelta);
    if constexpr (D == 2) {
      return meanCurvature2D(d);
    } else {
      return meanCurvature3D(d);
    }
  }

  template <class IteratorType>
  T gaussianCurvature(IteratorType &neighborIterator) {
    if constexpr (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }
    std::array<T, 9> d = bigStencilFromIterator(neighborIterator, gridDelta);
    return gaussianCurvature3D(d);
  }
  template <class IteratorType>
  std::array<T, 2> meanGaussianCurvature(IteratorType &neighborIterator) {
    if constexpr (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }
    std::array<T, 9> d = bigStencilFromIterator(neighborIterator, gridDelta);
    T mean = meanCurvature3D(d);
    T gauss = gaussianCurvature3D(d);
    return std::array<T, 2>{mean, gauss};
  }
  // TODO: Add Variation of normals.
};

} // namespace lsInternal

#endif // LS_CURVATURE_FORMULAS_HPP
