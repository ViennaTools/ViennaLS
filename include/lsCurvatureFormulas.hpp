#ifndef LS_CURVATURE_FORMULAS_HPP
#define LS_CURVATURE_FORMULAS_HPP

#include <lsDomain.hpp>

namespace lsInternal {

// Formulas for space curves and higher dimension can be found here:
// https://doi.org/10.1016/j.cagd.2005.06.005

// Calculates the Mean Curvature and/or the Gaussian Curvature of the level set
// function in 3D Calculates the Curvature of the level set function in 2D
// Requires an iterator that is big enough to calculate second order
// derivatives(e.g. hrleBoxIterator or hrleCartesianPlaneIterator)
//
//    Slice of a stencil: x axis horizontal y axis vertical
//
//        phi_np | phi_py | phi_pp
//        phi_nx | phi_0  | phi px
//        phi_nn | phi_ny | phi_nn

template <class T, int D> class baseCurvatureFormula {
public:
  template <class IteratorType>
  T meanCurvature(IteratorType &neighborIterator) {
    return lsDomain<T, D>::POS_VALUE;
  }

  template <class IteratorType>
  T gaussianCurvature(IteratorType &neighborIterator) {
    return lsDomain<T, D>::POS_VALUE;
  }

  template <class IteratorType>
  std::array<T, 2> meanGaussianCurvature(IteratorType &neighborIterator) {
    return std::array<T, 2>{lsDomain<T, D>::POS_VALUE,
                            lsDomain<T, D>::POS_VALUE};
  }
};

// Calculates the Curvature using the General Formula for implicit surfaces and
// standard first order approximations to the derivatives
template <class T, int D>
class curvaturGeneralFormula : public baseCurvatureFormula<T, D> {

private:
  T gridDelta;

  T twoGD = 0;
  T gDSq = 0;
  T fourGDsq = 0;

public:
  curvaturGeneralFormula(T mGD) : gridDelta(mGD) {
    twoGD = 1. / (2. * gridDelta);
    gDSq = 1. / (gridDelta * gridDelta);
    fourGDsq = 1. / (4. * gridDelta * gridDelta);
  }

  template <class IteratorType>
  T meanCurvature(IteratorType &neighborIterator) {

    // Stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      posUnit[i] = 1;
      negUnit[i] = -1;

      int second_pos = (i + 1) % D;

      T phi_0 = neighborIterator.getCenter().getValue() * gridDelta;

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = -1;
      negUnit[second_pos] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[i] = 0;
      negUnit[i] = 0;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = -1;

      // Calculate derivatives for Hessian
      d[i] = (phi_px - phi_nx) * twoGD;

      d[i + 3] = (phi_px - 2. * phi_0 + phi_nx) * gDSq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_pow3 = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_pow3 += d[i] * d[i];
    }
    norm_grad_pow3 = std::sqrt(norm_grad_pow3);

    norm_grad_pow3 = norm_grad_pow3 * norm_grad_pow3 * norm_grad_pow3;

    if (D == 2) {
      // Curvature formula for implicit planar Curve (2D)
      return (d[3] * d[1] * d[1] - 2. * d[1] * d[0] * d[6] +
              d[4] * d[0] * d[0]) /
             (norm_grad_pow3);
    }

    // Mean Curvature formula for implicit surfaces (3D)
    return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx +
        //    F_yy)
        (d[0] * d[0] * (d[4] + d[5]) + d[1] * d[1] * (d[3] + d[5]) +
         d[2] * d[2] * (d[3] + d[4]) +

         //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
         -2. * (d[0] * d[1] * d[6] + d[0] * d[2] * d[8] + d[1] * d[2] * d[7]))

        // /2*(F_x² + F_y² + F_z²)^(3/2)
        / (2. * norm_grad_pow3);
  }

  template <class IteratorType>
  T gaussianCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    if (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      posUnit[i] = 1;
      negUnit[i] = -1;

      int second_pos = (i + 1) % D;

      T phi_0 = neighborIterator.getCenter().getValue() * gridDelta;

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = -1;
      negUnit[second_pos] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[i] = 0;
      negUnit[i] = 0;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = -1;

      // Calculate derivatives for Hessian
      d[i] = (phi_px - phi_nx) * twoGD;

      d[i + 3] = (phi_px - 2. * phi_0 + phi_nx) * gDSq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_sqared = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_sqared += d[i] * d[i];
    }

    norm_grad_sqared = norm_grad_sqared * norm_grad_sqared;

    // Gaussian Curvature formula for implicit surfaces (3D)
    //  - (F_x²(F_yz²-F_yyF_zz)   +     F_y²(F_zx²-F_xxF_zz)     +
    //  F_z²(F_xy²-F_xxFyy) +
    return -(d[0] * d[0] * (d[7] * d[7] - d[4] * d[5]) +
             d[1] * d[1] * (d[8] * d[8] - d[3] * d[5]) +
             d[2] * d[2] * (d[6] * d[6] - d[3] * d[4]) +

             // 2*[F_xF_y(2F_zzF_xy-2F_zxF_yz) +
             2. * (d[0] * d[1] * (d[5] * d[6] - d[8] * d[7]) +

                   // F_xF_z(2F_yyF_zx-2F_xyF_yz) +
                   d[0] * d[2] * (d[4] * d[8] - d[6] * d[7]) +

                   // F_yF_z(2F_xxF_yz-2F_xyF_zx))
                   d[1] * d[2] * (d[3] * d[7] - d[6] * d[8])))

           // /(F_x² + F_y² + F_z²)²
           / norm_grad_sqared;
  }

  template <class IteratorType>
  std::array<T, 2> meanGaussianCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    if (D == 2) {
      lsMessage::getInstance()
          .addWarning(
              "2D structures do not have a Gaussian Curvature, use "
              "\"meanCurvature(IteratorType & neighborIterator)\" instead!")
          .print();
    }

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      posUnit[i] = 1;
      negUnit[i] = -1;

      int second_pos = (i + 1) % D;

      T phi_0 = neighborIterator.getCenter().getValue() * gridDelta;

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[second_pos] = -1;
      negUnit[second_pos] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue() * gridDelta;
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue() * gridDelta;

      posUnit[i] = 0;
      negUnit[i] = 0;

      posUnit[second_pos] = 1;
      negUnit[second_pos] = -1;

      // Calculate derivatives for Hessian
      d[i] = (phi_px - phi_nx) * twoGD;

      d[i + 3] = (phi_px - 2. * phi_0 + phi_nx) * gDSq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_pow3 = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_pow3 += d[i] * d[i];
    }

    T norm_grad_sqared = norm_grad_pow3 * norm_grad_pow3;

    norm_grad_pow3 = std::sqrt(norm_grad_pow3);

    norm_grad_pow3 = norm_grad_pow3 * norm_grad_pow3 * norm_grad_pow3;

    // Mean Curvature formula for implicit surfaces (3D)
    //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx +
    //    F_yy)
    T mean =
        (d[0] * d[0] * (d[4] + d[5]) + d[1] * d[1] * (d[3] + d[5]) +
         d[2] * d[2] * (d[3] + d[4]) +

         //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
         -2. * (d[0] * d[1] * d[6] + d[0] * d[2] * d[8] + d[1] * d[2] * d[7]))

        // /2*(F_x² + F_y² + F_z²)^(3/2)
        / (2. * norm_grad_pow3);
    // Gaussian Curvature formula for implicit surfaces (3D)
    //  - (F_x²(F_yz²-F_yyF_zz)   +     F_y²(F_zx²-F_xxF_zz)     +
    //  F_z²(F_xy²-F_xxFyy) +
    T gauss = -(d[0] * d[0] * (d[7] * d[7] - d[4] * d[5]) +
                d[1] * d[1] * (d[8] * d[8] - d[3] * d[5]) +
                d[2] * d[2] * (d[6] * d[6] - d[3] * d[4]) +

                // 2*[F_xF_y(2F_zzF_xy-2F_zxF_yz) +
                2. * (d[0] * d[1] * (d[5] * d[6] - d[8] * d[7]) +

                      // F_xF_z(2F_yyF_zx-2F_xyF_yz) +
                      d[0] * d[2] * (d[4] * d[8] - d[6] * d[7]) +

                      // F_yF_z(2F_xxF_yz-2F_xyF_zx))
                      d[1] * d[2] * (d[3] * d[7] - d[6] * d[8])))

              // /(F_x² + F_y² + F_z²)²
              / norm_grad_sqared;

    return std::array<T, 2>{mean, gauss};
  }
};

// Calculates the Curvature using the General Formula for implict surfaces and
// uses different first order approximations for D_x and D_xx derivatives
template <class T, int D>
class curvaturGeneralFormulaBigStencil : public baseCurvatureFormula<T, D> {

private:
  T gridDelta;

  T twoGD = 0;
  T fourGDsq = 0;
  T threeGDsq = 0;
  T fourGD = 0;

public:
  curvaturGeneralFormulaBigStencil(T mGD) : gridDelta(mGD) {

    twoGD = 1. / (2. * gridDelta);
    fourGD = 1. / (4. * gridDelta);
    threeGDsq = 1. / (3. * gridDelta * gridDelta);
    fourGDsq = 1. / (4. * gridDelta * gridDelta);
  }

  template <class IteratorType>
  T meanCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      int first_axis = i;
      int second_axis = (i + 1) % D;

      posUnit[first_axis] = 1;
      negUnit[first_axis] = -1;

      // get required ls values
      T phi_0 = neighborIterator.getCenter().getValue();

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = 1;
      negUnit[second_axis] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = -1;
      negUnit[second_axis] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[first_axis] = 0;
      negUnit[first_axis] = 0;

      posUnit[second_axis] = 1;
      negUnit[second_axis] = -1;

      T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

      d[i] = (phi_pp - phi_np + phi_pn - phi_nn) * fourGD;

      d[i + 3] = (phi_pp - 2. * phi_py + phi_np + phi_px - 2. * phi_0 + phi_nx +
                  phi_pn - 2. * phi_ny + phi_nn) *
                 threeGDsq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_pow3 = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_pow3 += d[i] * d[i];
    }
    norm_grad_pow3 = std::sqrt(norm_grad_pow3);

    norm_grad_pow3 = norm_grad_pow3 * norm_grad_pow3 * norm_grad_pow3;

    if (D == 2) {
      // Curvature formula for implicit planar Curve (2D)
      return (d[3] * d[1] * d[1] - 2. * d[1] * d[0] * d[6] +
              d[4] * d[0] * d[0]) /
             (norm_grad_pow3);
    }

    // Mean Curvature formula for implicit surfaces (3D)
    return
        //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx +
        //    F_yy)
        (d[0] * d[0] * (d[4] + d[5]) + d[1] * d[1] * (d[3] + d[5]) +
         d[2] * d[2] * (d[3] + d[4]) +

         //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
         -2. * (d[0] * d[1] * d[6] + d[0] * d[2] * d[8] + d[1] * d[2] * d[7]))

        // /2*(F_x² + F_y² + F_z²)^(3/2)
        / (2. * norm_grad_pow3);

    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
  }

  template <class IteratorType>
  T gaussianCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      int first_axis = i;
      int second_axis = (i + 1) % D;

      posUnit[first_axis] = 1;
      negUnit[first_axis] = -1;

      T phi_0 = neighborIterator.getCenter().getValue();

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = 1;
      negUnit[second_axis] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = -1;
      negUnit[second_axis] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[first_axis] = 0;
      negUnit[first_axis] = 0;

      posUnit[second_axis] = 1;
      negUnit[second_axis] = -1;

      T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

      d[i] = (phi_pp - phi_np + phi_pn - phi_nn) * fourGD;

      d[i + 3] = (phi_pp - 2. * phi_py + phi_np + phi_px - 2. * phi_0 + phi_nx +
                  phi_pn - 2. * phi_ny + phi_nn) *
                 threeGDsq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_sqared = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_sqared += d[i] * d[i];
    }

    norm_grad_sqared = norm_grad_sqared * norm_grad_sqared;

    // Gaussian Curvature formula for implicit surfaces (3D)
    //  - (F_x²(F_yz²-F_yyF_zz)   +     F_y²(F_zx²-F_xxF_zz)     +
    //  F_z²(F_xy²-F_xxFyy) +
    return -(d[0] * d[0] * (d[7] * d[7] - d[4] * d[5]) +
             d[1] * d[1] * (d[8] * d[8] - d[3] * d[5]) +
             d[2] * d[2] * (d[6] * d[6] - d[3] * d[4]) +

             // 2*[F_xF_y(2F_zzF_xy-2F_zxF_yz) +
             2. * (d[0] * d[1] * (d[5] * d[6] - d[8] * d[7]) +

                   // F_xF_z(2F_yyF_zx-2F_xyF_yz) +
                   d[0] * d[2] * (d[4] * d[8] - d[6] * d[7]) +

                   // F_yF_z(2F_xxF_yz-2F_xyF_zx))
                   d[1] * d[2] * (d[3] * d[7] - d[6] * d[8])))

           // /(F_x² + F_y² + F_z²)²
           / norm_grad_sqared;
  }
  template <class IteratorType>
  std::array<T, 2> meanGaussianCurvature(IteratorType &neighborIterator) {

    // stores the values of the discrete Hessian of the level set function
    //(F_x, F_y, F_z, F_xx, F_yy, F_zz, F_xy, F_yz, F_zx)
    std::array<T, 9> d;

    for (int i = 0; i < D; i++) {

      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);

      int first_axis = i;
      int second_axis = (i + 1) % D;

      posUnit[first_axis] = 1;
      negUnit[first_axis] = -1;

      T phi_0 = neighborIterator.getCenter().getValue();

      T phi_px = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nx = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = 1;
      negUnit[second_axis] = 1;

      T phi_pp = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_np = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[second_axis] = -1;
      negUnit[second_axis] = -1;

      T phi_pn = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_nn = neighborIterator.getNeighbor(negUnit).getValue();

      posUnit[first_axis] = 0;
      negUnit[first_axis] = 0;

      posUnit[second_axis] = 1;
      negUnit[second_axis] = -1;

      T phi_py = neighborIterator.getNeighbor(posUnit).getValue();
      T phi_ny = neighborIterator.getNeighbor(negUnit).getValue();

      d[i] = (phi_pp - phi_np + phi_pn - phi_nn) * fourGD;

      d[i + 3] = (phi_pp - 2. * phi_py + phi_np + phi_px - 2. * phi_0 + phi_nx +
                  phi_pn - 2. * phi_ny + phi_nn) *
                 threeGDsq;

      d[i + 6] = (phi_pp - phi_pn - phi_np + phi_nn) * fourGDsq;
    }

    T norm_grad_pow3 = 0.;
    for (int i = 0; i < D; i++) {

      norm_grad_pow3 += d[i] * d[i];
    }

    T norm_grad_sqared = norm_grad_pow3 * norm_grad_pow3;

    norm_grad_pow3 = std::sqrt(norm_grad_pow3);

    norm_grad_pow3 = norm_grad_pow3 * norm_grad_pow3 * norm_grad_pow3;

    // Mean Curvature formula for implicit surfaces (3D)
    //    F_x²(f_yy + F_zz)  +    F_y²(F_xx + F_zz)    +     F_z²(F_xx +
    //    F_yy)
    T mean =
        (d[0] * d[0] * (d[4] + d[5]) + d[1] * d[1] * (d[3] + d[5]) +
         d[2] * d[2] * (d[3] + d[4]) +

         //-2*[F_xF_yF_xy   +   F_xF_zF_xz   +   F_yF_zF_yz]
         -2. * (d[0] * d[1] * d[6] + d[0] * d[2] * d[8] + d[1] * d[2] * d[7]))

        // /2*(F_x² + F_y² + F_z²)^(3/2)
        / (2. * norm_grad_pow3);
    // Gaussian Curvature formula for implicit surfaces (3D)
    //  - (F_x²(F_yz²-F_yyF_zz)   +     F_y²(F_zx²-F_xxF_zz)     +
    //  F_z²(F_xy²-F_xxFyy) +
    T gauss = -(d[0] * d[0] * (d[7] * d[7] - d[4] * d[5]) +
                d[1] * d[1] * (d[8] * d[8] - d[3] * d[5]) +
                d[2] * d[2] * (d[6] * d[6] - d[3] * d[4]) +

                // 2*[F_xF_y(2F_zzF_xy-2F_zxF_yz) +
                2. * (d[0] * d[1] * (d[5] * d[6] - d[8] * d[7]) +

                      // F_xF_z(2F_yyF_zx-2F_xyF_yz) +
                      d[0] * d[2] * (d[4] * d[8] - d[6] * d[7]) +

                      // F_yF_z(2F_xxF_yz-2F_xyF_zx))
                      d[1] * d[2] * (d[3] * d[7] - d[6] * d[8])))

              // /(F_x² + F_y² + F_z²)²
              / norm_grad_sqared;

    return std::array<T, 2>{mean, gauss};
  }

  // TODO: Add Variation of normals.
};

} // namespace lsInternal

#endif // LS_CURVATURE_FORMULAS_HPP
