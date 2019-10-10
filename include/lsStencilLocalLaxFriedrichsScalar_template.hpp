#ifndef LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_TEMPLATE_HPP
#define LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_TEMPLATE_HPP

#include <hrleSquareIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain_template.hpp>
#include <lsExpand_template.hpp>
#include <lsFiniteDifferences_template.hpp>

namespace lsInternal {

/// Stencil Local Lax Friedrichs Integration Scheme.
/// It uses a stencil of order around active points, in order to
/// evaluate dissipation values for each point, taking into account
/// the mathematical nature of the speed function.
/// see Toifl et al., 2019. ISBN: 978-1-7281-0938-1
template <class T, int D, int order> class lsStencilLocalLaxFriedrichsScalar {
  lsDomain<T, D> &levelSet;
  const DifferentiationSchemeEnum finiteDifferenceScheme =
      DifferentiationSchemeEnum::FIRST_ORDER;
  // const int stencilOrder = 1;
  hrleSquareIterator<hrleDomain<T, D>> neighborIterator;

  // Final dissipation coefficients that are used by the time integrator. If
  // D==2 last entries are 0.
  hrleVectorType<T, 3> finalAlphas;
  hrleVectorType<T, 3> gridDeltas;
  const unsigned numStencilPoints;

  static T pow2(const T &value) { return value * value; }

  hrleVectorType<T, D>
  calculateNormal(const hrleVectorType<hrleIndexType, D> &offset) {
    hrleVectorType<T, D> normal;
    const int startIndex = -1;
    T modulus = 0.;

    for (unsigned i = 0; i < D; ++i) {
      hrleVectorType<hrleIndexType, D> index(0);
      std::vector<T> values;
      for (unsigned j = 0; j < 3; ++j) {
        index[i] = startIndex + j;
        values.push_back(
            neighborIterator.getNeighbor(offset + index).getValue());
      }
      normal[i] = lsFiniteDifferences<T>::calculateGradient(
          &(values[0]), levelSet.getGrid().getGridDelta());
      modulus += normal[i] * normal[i];
    }
    modulus = std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      normal[i] /= modulus;
    }
    return normal;
  }

  hrleVectorType<T, D>
  calculateGradient(const hrleVectorType<hrleIndexType, D> &offset) {
    hrleVectorType<T, D> gradient;

    const unsigned numValues =
        static_cast<unsigned>(finiteDifferenceScheme) + 2;
    const int startIndex = -std::floor(numValues / 2);

    for (unsigned i = 0; i < D; ++i) {
      hrleVectorType<hrleIndexType, D> index(hrleIndexType(0));
      std::vector<T> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values.push_back(
            neighborIterator.getNeighbor(offset + index).getValue());
      }

      if (finiteDifferenceScheme == DifferentiationSchemeEnum::FIRST_ORDER) {
        gradient[i] =
            lsFiniteDifferences<T, DifferentiationSchemeEnum::FIRST_ORDER>::
                calculateGradient(&(values[0]),
                                  levelSet.getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO3) {
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO3>::
            calculateGradient(&(values[0]), levelSet.getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO5)
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO5>::
            calculateGradient(&(values[0]), levelSet.getGrid().getGridDelta());
    }

    return gradient;
  }

  hrleVectorType<T, D> calculateGradientDiff() {
    hrleVectorType<T, D> gradient;

    const unsigned numValues =
        static_cast<unsigned>(finiteDifferenceScheme) + 2;
    const int startIndex = -std::floor(numValues / 2);

    for (unsigned i = 0; i < D; ++i) {
      hrleVectorType<hrleIndexType, D> index(hrleIndexType(0));
      std::vector<T> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }

      if (finiteDifferenceScheme == DifferentiationSchemeEnum::FIRST_ORDER) {
        gradient[i] =
            lsFiniteDifferences<T, DifferentiationSchemeEnum::FIRST_ORDER>::
                calculateGradientDiff(&(values[0]),
                                      levelSet.getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO3) {
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO3>::
            calculateGradientDiff(&(values[0]),
                                  levelSet.getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO5)
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO5>::
            calculateGradientDiff(&(values[0]),
                                  levelSet.getGrid().getGridDelta());
    }

    return gradient;
  }

public:
  const hrleVectorType<T, 3> &getFinalAlphas() const { return finalAlphas; }
  const hrleVectorType<T, 3> &getDeltas() const { return gridDeltas; }

  static void prepareLS(lsDomain<T, D> &passedlsDomain) {
    // Expansion of sparse field must depend on spatial derivative order
    // AND  slf stencil order! --> currently assume scheme = 3rd order always
    lsExpand<T, D>(passedlsDomain).apply(2 * order + 4);
  }

  lsStencilLocalLaxFriedrichsScalar(
      lsDomain<T, D> &passedlsDomain,
      DifferentiationSchemeEnum scheme = DifferentiationSchemeEnum::FIRST_ORDER)
      : levelSet(passedlsDomain), finiteDifferenceScheme(scheme),
        neighborIterator(hrleSquareIterator<hrleDomain<T, D>>(
            levelSet.getDomain(), static_cast<unsigned>(scheme) + 1 + order)),
        numStencilPoints(std::pow(2 * order + 1, D)) {
    for (int i = 0; i < 3; ++i) {
      finalAlphas[i] = 0;
      gridDeltas[i] = 0;
    }
  }

  T operator()(const hrleVectorType<hrleIndexType, D> &indices,
               lsVelocityField<T> *velocities, int material) {
    auto &grid = levelSet.getGrid();
    double gridDelta = grid.getGridDelta();

    hrleVectorType<double, 3> coordinate(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    hrleVectorType<double, 3> vectorVelocity =
        velocities->getVectorVelocity(coordinate, material);
    double scalarVelocity = velocities->getScalarVelocity(coordinate, material);

    // if there is a vector velocity, we need to project it onto a scalar
    // velocity first using its normal vector
    if (vectorVelocity != hrleVectorType<T, D>(0., 0., 0.)) {
      hrleVectorType<T, D> n;
      T denominator = 0; // normal modulus
      for (unsigned i = 0; i < D; i++) {
        int index[3] = {i == 0, i == 1, i == 2}; // unit vector in direction
        // normal vector calculation
        T pos = neighborIterator.getNeighbor(index).getValue() -
                neighborIterator.getCenter().getValue();
        index[i] = -index[i];
        T neg = neighborIterator.getCenter().getValue() -
                neighborIterator.getNeighbor(index).getValue();
        n[i] = (pos + neg) * 0.5;
        // normalise normal vector
        denominator += n[i] * n[i];
      }
      denominator = std::sqrt(denominator);
      // now calculate scalar product of normal vector with velocity
      for (unsigned i = 0; i < D; ++i) {
        scalarVelocity += vectorVelocity[i] * n[i] / denominator;
      }
    }

    if (scalarVelocity == T(0)) {
      return 0;
    }

    // const bool DEBUG = false;

    // typename LevelSetType::neighbor_stencil n(LS, it, stencilOrder);
    // std::vector<typename LevelSetType::star_stencil> stars =
    //     n.star_stencils(LevelSetType::differentiation_scheme::FIRST_ORDER);
    // typename LevelSetType::star_stencil center = stars[n.get_center_index()];

    // hrleVectorType<T, D> dxs = center.getDx();
    for (unsigned i = 0; i < D; ++i) {
      gridDeltas[i] = gridDelta;
    }

    // if (DEBUG) {
    //   for (auto star : stars) {
    //     std::cout << star.position() << std::endl;
    //     std::cout << star.center().value() << std::endl;
    //     std::cout << star.center().active_pt_id() << std::endl;
    //   }
    //   std::cout << std::endl;
    // }

    T hamiltonian =
        NormL2(calculateGradient(hrleVectorType<hrleIndexType, D>(0))) *
        scalarVelocity;
    T dissipation = 0.; // dissipation

    // if (DEBUG)
    //   std::cout << "H = " << hamiltonian << std::endl;

    // dissipation block
    {
      // reserve alphas for all points in local stencil
      std::vector<hrleVectorType<T, D>> alphas;
      alphas.reserve(numStencilPoints);

      hrleVectorType<hrleIndexType, D> currentIndex(-order);
      for (size_t i = 0; i < numStencilPoints; ++i) {
        hrleVectorType<T, D> alpha;
        hrleVectorType<T, 3> normal(calculateNormal(currentIndex));
        if (D == 2)
          normal[2] = 0;

        // Check for corrupted normal
        if ((std::abs(normal[0]) < 1e-6) && (std::abs(normal[1]) < 1e-6) &&
            (std::abs(normal[2]) < 1e-6)) {
          alphas.push_back(hrleVectorType<T, D>(0.0));
          continue;
        }

        hrleVectorType<T, 3> normal_n = normal;
        hrleVectorType<T, 3> normal_p = normal;

        hrleVectorType<T, D> velocityDelta(T(0));
        const T DN = 1e-4;

        for (int k = 0; k < D; ++k) {

          normal_p[k] -= DN; // p=previous
          normal_n[k] += DN; // n==next

          T vp = velocities->getScalarVelocity(coordinate, material, normal_p);
          T vn = velocities->getScalarVelocity(coordinate, material, normal_n);
          // central difference
          velocityDelta[k] = (vn - vp) / (2.0 * DN);

          normal_p[k] += DN;
          normal_n[k] -= DN;
        }

        // determine \partial H / \partial phi_l
        for (int k = 0; k < D; ++k) { // iterate over dimensions

          // Monti term
          T monti = 0;
          // Toifl Quell term
          T toifl = 0;

          hrleVectorType<T, D> gradient = calculateGradient(currentIndex);

          for (int j = 0; j < D - 1; ++j) { // phi_p**2 + phi_q**2
            int idx = (k + 1 + j) % D;
            monti += gradient[idx] * gradient[idx];
            toifl += gradient[idx] * velocityDelta[idx];
          }
          // denominator: |grad(phi)|^2
          T denominator = Norm2(gradient);
          monti *= velocityDelta[k] / denominator;
          toifl *= -gradient[k] / denominator;

          // Osher (constant V) term
          T osher = velocities->getScalarVelocity(coordinate, material, normal);

          // Total derivative is sum of terms given above
          alpha[k] = std::fabs(monti + toifl + osher);
        }

        alphas.push_back(alpha);

        // increment current index
        {
          int dim = 0;
          for (; dim < D - 1; ++dim) {
            if (currentIndex[dim] < order)
              break;
            currentIndex[dim] = -order;
          }
          ++currentIndex[dim];
        }
      }

      // determine max alphas for every axis
      hrleVectorType<T, D> maxal;
      hrleVectorType<T, D> gradientDiff = calculateGradientDiff();
      for (int d = 0; d < D; ++d) {
        maxal[d] = 0;
        std::vector<T> alpha_comp;

        for (size_t i = 0; i < numStencilPoints; ++i) {
          alpha_comp.push_back(alphas[i][d]);
        }
        maxal[d] = *std::max_element(alpha_comp.begin(), alpha_comp.end());

        finalAlphas[d] = maxal[d];

        dissipation += maxal[d] * gradientDiff[d];
      }
      //   if (DEBUG)
      //     std::cout << "max(alpha) = " << maxal << std::endl << std::endl;
      //   if (DEBUG)
      //     std::cout << "D = " << dissipation << std::endl;
      //
      //   if (DEBUG)
      //     std::cout << "diff = " << center.gradient_diff() << "\n";
      // }
      //
      // if (DEBUG)
      //   std::cout << "H-D = " << hamiltonian - dissipation << std::endl;

      return hamiltonian - dissipation;
    }
  }
};

} // namespace lsInternal

#endif // LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_TEMPLATE_HPP
