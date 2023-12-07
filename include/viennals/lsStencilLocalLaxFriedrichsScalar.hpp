#ifndef LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_HPP
#define LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_HPP

#include <hrleSparseBoxIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>

namespace lsInternal {

/// Stencil Local Lax Friedrichs Integration Scheme.
/// It uses a stencil of order around active points, in order to
/// evaluate dissipation values for each point, taking into account
/// the mathematical nature of the speed function.
/// see Toifl et al., 2019. ISBN: 978-1-7281-0938-1;
/// DOI: 10.1109/SISPAD.2019.8870443
template <class T, int D, int order> class lsStencilLocalLaxFriedrichsScalar {
  using LevelSetType = lsSmartPointer<lsDomain<T, D>>;
  using LevelSetsType = std::vector<LevelSetType>;

  LevelSetType levelSet;
  lsSmartPointer<lsVelocityField<T>> velocities;
  const DifferentiationSchemeEnum finiteDifferenceScheme =
      DifferentiationSchemeEnum::FIRST_ORDER;
  hrleSparseBoxIterator<hrleDomain<T, D>> neighborIterator;
  const double alphaFactor;
  const double normalEpsilon =
      std::cbrt(std::numeric_limits<double>::epsilon());

  // Final dissipation coefficients that are used by the time integrator. If
  // D==2 last entries are 0.
  hrleVectorType<T, 3> finalAlphas;
  const unsigned numStencilPoints;

  static T pow2(const T &value) { return value * value; }

  hrleVectorType<T, D>
  calculateNormal(const hrleVectorType<hrleIndexType, D> &offset) {
    hrleVectorType<T, D> normal;
    const int startIndex = -1;
    T modulus = 0.;

    for (unsigned i = 0; i < D; ++i) {
      hrleVectorType<hrleIndexType, D> index(offset);
      std::vector<T> values;
      for (unsigned j = 0; j < 3; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }
      normal[i] = lsFiniteDifferences<T>::calculateGradient(
          &(values[0]), levelSet->getGrid().getGridDelta());
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
        lsFiniteDifferences<T>::getNumberOfValues(finiteDifferenceScheme);
    const int startIndex = -std::floor(numValues / 2);

    for (unsigned i = 0; i < D; ++i) {
      hrleVectorType<hrleIndexType, D> index(offset);
      std::vector<T> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }

      if (finiteDifferenceScheme == DifferentiationSchemeEnum::FIRST_ORDER) {
        gradient[i] =
            lsFiniteDifferences<T, DifferentiationSchemeEnum::FIRST_ORDER>::
                calculateGradient(&(values[0]),
                                  levelSet->getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO3) {
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO3>::
            calculateGradient(&(values[0]), levelSet->getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO5)
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO5>::
            calculateGradient(&(values[0]), levelSet->getGrid().getGridDelta());
    }

    return gradient;
  }

  hrleVectorType<T, D> calculateGradientDiff() {
    hrleVectorType<T, D> gradient;

    const unsigned numValues =
        lsFiniteDifferences<T>::getNumberOfValues(finiteDifferenceScheme);
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
                                      levelSet->getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO3) {
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO3>::
            calculateGradientDiff(&(values[0]),
                                  levelSet->getGrid().getGridDelta());
      } else if (finiteDifferenceScheme == DifferentiationSchemeEnum::WENO5)
        gradient[i] = lsFiniteDifferences<T, DifferentiationSchemeEnum::WENO5>::
            calculateGradientDiff(&(values[0]),
                                  levelSet->getGrid().getGridDelta());
    }

    return gradient;
  }

public:
  const hrleVectorType<T, 3> &getFinalAlphas() const { return finalAlphas; }

  static void prepareLS(LevelSetType passedlsDomain) {
    // Expansion of sparse field must depend on spatial derivative order
    // AND  slf stencil order! --> currently assume scheme = 3rd order always
    lsExpand<T, D>(passedlsDomain, 2 * (order + 1) + 4).apply();
  }

  lsStencilLocalLaxFriedrichsScalar(
      LevelSetType passedlsDomain, lsSmartPointer<lsVelocityField<T>> vel,
      double a = 1.0,
      DifferentiationSchemeEnum scheme = DifferentiationSchemeEnum::FIRST_ORDER)
      : levelSet(passedlsDomain), velocities(vel),
        finiteDifferenceScheme(scheme),
        neighborIterator(hrleSparseBoxIterator<hrleDomain<T, D>>(
            levelSet->getDomain(), static_cast<unsigned>(scheme) + 1 + order)),
        alphaFactor(a), numStencilPoints(std::pow(2 * order + 1, D)) {

    for (int i = 0; i < 3; ++i) {
      finalAlphas[i] = 0;
    }
  }

  T operator()(const hrleVectorType<hrleIndexType, D> &indices, int material) {
    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    hrleVectorType<T, 3> coordinate(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    // convert coordinate to std array for interface
    std::array<T, 3> coordArray = {coordinate[0], coordinate[1], coordinate[2]};

    // if there is a vector velocity, we need to project it onto a scalar
    // velocity first using its normal vector
    // /*if (vectorVelocity != std::array<T, 3>({}))*/ {
    std::array<T, 3> normalVector = {};
    T denominator = 0; // normal modulus
    for (unsigned i = 0; i < D; i++) {
      hrleVectorType<T, 3> neighborIndex(T(0));
      neighborIndex[i] = 1;
      // normal vector calculation
      T pos = neighborIterator.getNeighbor(neighborIndex).getValue() -
              neighborIterator.getCenter().getValue();
      T neg = neighborIterator.getCenter().getValue() -
              neighborIterator.getNeighbor(-neighborIndex).getValue();
      normalVector[i] = (pos + neg) * 0.5;
      // normalise normal vector
      denominator += normalVector[i] * normalVector[i];
    }
    denominator = std::sqrt(denominator);

    for (unsigned i = 0; i < D; ++i) {
      normalVector[i] /= denominator;
    }

    double scalarVelocity = velocities->getScalarVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());
    std::array<T, 3> vectorVelocity = velocities->getVectorVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());

    // now calculate scalar product of normal vector with velocity
    for (unsigned i = 0; i < D; ++i) {
      scalarVelocity += vectorVelocity[i] * normalVector[i];
    }

    if (scalarVelocity == T(0)) {
      return 0;
    }

    T hamiltonian =
        NormL2(calculateGradient(hrleVectorType<hrleIndexType, D>(0))) *
        scalarVelocity;
    T dissipation = 0.; // dissipation

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
          alphas.push_back(hrleVectorType<T, D>(T(0)));
          continue;
        }

        std::array<T, 3> normal_p = {normal[0], normal[1], normal[2]};
        std::array<T, 3> normal_n = {normal[0], normal[1], normal[2]};

        hrleVectorType<T, D> velocityDelta(T(0));

        // get local velocity
        std::array<T, 3> localCoordArray = coordArray;
        for (unsigned dir = 0; dir < D; ++dir)
          localCoordArray[dir] += currentIndex[dir];

        T localScalarVelocity = velocities->getScalarVelocity(
            localCoordArray, material, normal_p,
            neighborIterator.getCenter().getPointId());
        std::array<T, 3> localVectorVelocity = velocities->getVectorVelocity(
            localCoordArray, material, normal_p,
            neighborIterator.getCenter().getPointId());
        // now calculate scalar product of normal vector with velocity
        for (unsigned i = 0; i < D; ++i) {
          localScalarVelocity += localVectorVelocity[i] * normal[i];
        }

        const T DN = std::abs(normalEpsilon * scalarVelocity);

        for (int k = 0; k < D; ++k) {

          normal_p[k] -= DN; // p=previous
          normal_n[k] += DN; // n==next

          T vp = velocities->getScalarVelocity(
              localCoordArray, material, normal_p,
              neighborIterator.getCenter().getPointId());
          T vn = velocities->getScalarVelocity(
              localCoordArray, material, normal_n,
              neighborIterator.getCenter().getPointId());
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
          T osher = localScalarVelocity * normal[k];
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
      hrleVectorType<T, D> gradientDiff = calculateGradientDiff();
      for (int d = 0; d < D; ++d) {
        T maxAlpha = 0;

        for (size_t i = 0; i < numStencilPoints; ++i) {
          maxAlpha = std::max(maxAlpha, alphas[i][d]);
        }

        finalAlphas[d] = maxAlpha;
        dissipation += maxAlpha * gradientDiff[d];
      }

      return hamiltonian - dissipation;
    }
  }
};

namespace advect {
template <
    class IntegrationSchemeType, class T, int D,
    lsConcepts::IsSame<IntegrationSchemeType,
                       lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>> =
        lsConcepts::assignable>
void reduceTimeStepHamiltonJacobi(IntegrationSchemeType &scheme,
                                  double &MaxTimeStep,
                                  hrleCoordType gridDelta) {
  const double alpha_maxCFL = 1.0;
  // second time step test, based on alphas
  hrleVectorType<T, 3> alphas = scheme.getFinalAlphas();

  double timeStep = 0;
  for (int i = 0; i < D; ++i) {
    timeStep += alphas[i] / gridDelta;
  }

  timeStep = alpha_maxCFL / timeStep;
  MaxTimeStep = std::min(timeStep, MaxTimeStep);
}
} // namespace advect
} // namespace lsInternal

/// This function creates the specialized layer wrapping which
/// produces better results for the SSLF integration scheme.
/// isDepo must contain whether the corresponding level sets
/// are used for deposition or not.
/// This function assumes that the layers where deposition is
/// not possible represent only masking layers and are not
/// wrapped around other depo layers.
/// Hint: Sometimes it is useful to introduce a new
/// mask layer combining all masking materials and remove it
/// after advection instead of trying to deal with numerous
/// separate masking layers.
template <class T, int D>
void lsPrepareStencilLocalLaxFriedrichs(
    std::vector<lsSmartPointer<lsDomain<T, D>>> &levelSets,
    std::vector<bool> isDepo) {
  if (isDepo.size() < levelSets.size()) {
    lsMessage::getInstance()
        .addWarning(
            "lsPrepareStencilLocalLaxFriedrichs: isDepo does not have enough "
            "elements. Assuming all higher layers are not depo layers.")
        .print();
  }
  // always resize, so it has the correct number of elements
  isDepo.resize(levelSets.size(), false);

  // Begin with biggest level set (top LS wrapped around all others)
  auto layerIt = levelSets.rbegin();

  // make empty LS which will contain the final top layer
  auto finalTop = lsSmartPointer<lsDomain<T, D>>::New(levelSets[0]->getGrid());

  bool layerAboveIsDepo = false;
  auto depoIt = isDepo.rbegin();

  // subtract each first non-depo layer below a depo layer from the latter
  // then union these results to the final layer
  for (auto maskIt = levelSets.rbegin(); maskIt != levelSets.rend(); ++maskIt) {
    if (isDepo[*depoIt]) {
      // this layer will be deposited on
      if (!layerAboveIsDepo) {
        layerIt = maskIt;
        layerAboveIsDepo = true;
      }
    } else {
      // no depo on this layer
      if (layerAboveIsDepo) {
        auto layerAbove = lsSmartPointer<lsDomain<T, D>>::New(*layerIt);
        lsBooleanOperation<T, D>(layerAbove, *maskIt,
                                 lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
            .apply();

        lsBooleanOperation<T, D>(finalTop, layerAbove,
                                 lsBooleanOperationEnum::UNION)
            .apply();
        lsPrune<T, D>(finalTop).apply();
      }
      layerAboveIsDepo = false;
    }
    ++depoIt;
  }

  // If the lowest layer is a depo substrate, add it to the final top layer
  if (layerAboveIsDepo) {
    lsBooleanOperation<T, D>(finalTop, *layerIt, lsBooleanOperationEnum::UNION)
        .apply();
  }

  lsBooleanOperation<T, D>(finalTop, lsBooleanOperationEnum::INVERT).apply();
  levelSets.back()->deepCopy(finalTop);
}

/// After advection using the SLLF layer wrapping approach is done,
/// restore the original layer wrapping used everywhere else.
template <class T, int D>
void lsFinalizeStencilLocalLaxFriedrichs(
    std::vector<lsSmartPointer<lsDomain<T, D>>> &levelSets) {
  auto layerIt = levelSets.rbegin();
  auto lastIt = ++levelSets.rbegin();

  lsBooleanOperation<T, D>(*layerIt, lsBooleanOperationEnum::INVERT).apply();
  lsBooleanOperation<T, D>(*layerIt, *lastIt, lsBooleanOperationEnum::UNION)
      .apply();
}

#endif // LS_STENCIL_LOCAL_LACHS_FRIEDRICHS_SCALAR_HPP
