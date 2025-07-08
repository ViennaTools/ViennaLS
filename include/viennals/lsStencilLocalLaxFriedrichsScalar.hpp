#pragma once

#include <hrleSparseBoxIterator.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>
#include <lsVelocityField.hpp>

#include <vcVectorType.hpp>

namespace lsInternal {

using namespace viennacore;

/// Stencil Local Lax Friedrichs Integration Scheme.
/// It uses a stencil of order around active points, in order to
/// evaluate dissipation values for each point, taking into account
/// the mathematical nature of the speed function.
/// see Toifl et al., 2019. ISBN: 978-1-7281-0938-1;
/// DOI: 10.1109/SISPAD.2019.8870443
template <class T, int D, int order,
          DifferentiationSchemeEnum finiteDifferenceScheme =
              DifferentiationSchemeEnum::FIRST_ORDER>
class StencilLocalLaxFriedrichsScalar {
  using LevelSetType = SmartPointer<viennals::Domain<T, D>>;
  using LevelSetsType = std::vector<LevelSetType>;

  LevelSetType levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;
  viennahrle::SparseBoxIterator<viennahrle::Domain<T, D>,
                                static_cast<int>(finiteDifferenceScheme) + 1 +
                                    order>
      neighborIterator;
  const double alphaFactor;
  const double normalEpsilon =
      std::cbrt(std::numeric_limits<double>::epsilon());

  // Final dissipation coefficients that are used by the time integrator. If
  // D==2 last entries are 0.
  Vec3D<T> finalAlphas;
  static constexpr unsigned numStencilPoints = hrleUtil::pow(2 * order + 1, D);
  static double maxDissipation; // default: std::numeric_limits<double>::max();

  static T pow2(const T &value) { return value * value; }

  Vec3D<T> calculateNormal(const viennahrle::Index<D> &offset) {
    Vec3D<T> normal = {0.0, 0.0, 0.0};
    constexpr int startIndex = -1;
    T modulus = 0.;

    for (unsigned i = 0; i < D; ++i) {
      viennahrle::Index<D> index(offset);
      std::vector<T> values;
      for (unsigned j = 0; j < 3; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }
      normal[i] = FiniteDifferences<T>::calculateGradient(
          values.data(), levelSet->getGrid().getGridDelta());
      modulus += normal[i] * normal[i];
    }
    modulus = std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      normal[i] /= modulus;
    }
    return normal;
  }

  VectorType<T, D> calculateGradient(const viennahrle::Index<D> &offset) {
    VectorType<T, D> gradient;

    constexpr unsigned numValues =
        FiniteDifferences<T, finiteDifferenceScheme>::getNumberOfValues();
    const int startIndex = -std::floor(numValues / 2);

    for (unsigned i = 0; i < D; ++i) {
      viennahrle::Index<D> index(offset);
      std::vector<T> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }

      gradient[i] =
          FiniteDifferences<T, finiteDifferenceScheme>::calculateGradient(
              values.data(), levelSet->getGrid().getGridDelta());
    }

    return gradient;
  }

  VectorType<T, D> calculateGradientDiff() {
    VectorType<T, D> gradient;

    const unsigned numValues =
        FiniteDifferences<T, finiteDifferenceScheme>::getNumberOfValues();
    const int startIndex = -std::floor(numValues / 2);

    for (unsigned i = 0; i < D; ++i) {
      viennahrle::Index<D> index(0);
      std::vector<T> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values.push_back(neighborIterator.getNeighbor(index).getValue());
      }

      gradient[i] =
          FiniteDifferences<T, finiteDifferenceScheme>::calculateGradientDiff(
              &(values[0]), levelSet->getGrid().getGridDelta());
    }

    return gradient;
  }

public:
  static void prepareLS(LevelSetType passedlsDomain) {
    // Expansion of sparse field must depend on spatial derivative order
    // AND  slf stencil order! --> currently assume scheme = 3rd order always
    viennals::Expand<T, D>(passedlsDomain, 2 * (order + 1) + 4).apply();
  }

  StencilLocalLaxFriedrichsScalar(LevelSetType passedlsDomain,
                                  SmartPointer<viennals::VelocityField<T>> vel,
                                  double a = 1.0)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain()), alphaFactor(a) {
    for (int i = 0; i < 3; ++i) {
      finalAlphas[i] = 0;
    }
  }

  static void setMaxDissipation(double maxDiss) { maxDissipation = maxDiss; }

  std::pair<T, T> operator()(const viennahrle::Index<D> &indices,
                             int material) {
    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    Vec3D<T> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    // convert coordinate to std array for interface
    Vec3D<T> coordArray{coordinate[0], coordinate[1], coordinate[2]};

    // if there is a vector velocity, we need to project it onto a scalar
    // velocity first using its normal vector
    // /*if (vectorVelocity != Vec3D<T>({}))*/ {
    Vec3D<T> normalVector;
    T denominator = 0; // normal modulus
    for (unsigned i = 0; i < D; i++) {
      viennahrle::Index<D> neighborIndex(0);
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
    auto vectorVelocity = velocities->getVectorVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());

    // now calculate scalar product of normal vector with velocity
    for (unsigned i = 0; i < D; ++i) {
      scalarVelocity += vectorVelocity[i] * normalVector[i];
    }

    if (scalarVelocity == T(0)) {
      return {0, 0};
    }

    T hamiltonian =
        Norm(calculateGradient(viennahrle::Index<D>(0))) * scalarVelocity;
    T dissipation = 0.; // dissipation

    // dissipation block
    {
      // reserve alphas for all points in local stencil
      std::vector<VectorType<T, D>> alphas;
      alphas.reserve(numStencilPoints);

      viennahrle::Index<D> currentIndex(-order);
      for (size_t i = 0; i < numStencilPoints; ++i) {
        VectorType<T, D> alpha;
        Vec3D<T> normal = calculateNormal(currentIndex);

        // Check for corrupted normal
        if ((std::abs(normal[0]) < 1e-6) && (std::abs(normal[1]) < 1e-6) &&
            (std::abs(normal[2]) < 1e-6)) {
          continue;
        }

        Vec3D<T> normal_p{normal[0], normal[1], normal[2]};
        Vec3D<T> normal_n{normal[0], normal[1], normal[2]};

        VectorType<T, D> velocityDelta;
        std::fill(velocityDelta.begin(), velocityDelta.end(), 0.);

        // get local velocity
        Vec3D<T> localCoordArray = coordArray;
        for (unsigned dir = 0; dir < D; ++dir)
          localCoordArray[dir] += currentIndex[dir];

        T localScalarVelocity = velocities->getScalarVelocity(
            localCoordArray, material, normal_p,
            neighborIterator.getCenter().getPointId());
        Vec3D<T> localVectorVelocity = velocities->getVectorVelocity(
            localCoordArray, material, normal_p,
            neighborIterator.getCenter().getPointId());
        // now calculate scalar product of normal vector with velocity
        for (unsigned dir = 0; dir < D; ++dir) {
          localScalarVelocity += localVectorVelocity[dir] * normal[dir];
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

          VectorType<T, D> gradient = calculateGradient(currentIndex);

          for (int j = 0; j < D - 1; ++j) { // phi_p**2 + phi_q**2
            int idx = (k + 1 + j) % D;
            monti += gradient[idx] * gradient[idx];
            toifl += gradient[idx] * velocityDelta[idx];
          }
          // denominator: |grad(phi)|^2
          T denom = DotProduct(gradient, gradient);
          monti *= velocityDelta[k] / denom;
          toifl *= -gradient[k] / denom;

          // Osher (constant V) term
          T osher = localScalarVelocity * normal[k];
          // Total derivative is sum of terms given above
          alpha[k] = std::fabs(monti + toifl + osher);

          if (alpha[k] > maxDissipation) {
            alpha[k] = 0.;
          }
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
      VectorType<T, D> gradientDiff = calculateGradientDiff();
      for (int d = 0; d < D; ++d) {
        T maxAlpha = 0;

        for (size_t i = 0; i < numStencilPoints; ++i) {
          maxAlpha = std::max(maxAlpha, alphas[i][d]);
        }

        finalAlphas[d] = maxAlpha;
        dissipation += maxAlpha * gradientDiff[d];
      }

      return {hamiltonian, dissipation};
    }
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep,
                                    double gridDelta) const {
    constexpr double alpha_maxCFL = 1.0;
    // second time step test, based on alphas

    double timeStep = 0;
    for (int i = 0; i < D; ++i) {
      timeStep += finalAlphas[i] / gridDelta;
    }

    timeStep = alpha_maxCFL / timeStep;
    MaxTimeStep = std::min(timeStep, MaxTimeStep);
  }
};

template <class T, int D, int order,
          DifferentiationSchemeEnum finiteDifferenceScheme>
double StencilLocalLaxFriedrichsScalar<T, D, order,
                                       finiteDifferenceScheme>::maxDissipation =
    std::numeric_limits<double>::max();

} // namespace lsInternal

namespace viennals {

using namespace viennacore;

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
void PrepareStencilLocalLaxFriedrichs(
    std::vector<SmartPointer<Domain<T, D>>> &levelSets,
    std::vector<bool> isDepo) {
  if (isDepo.size() < levelSets.size()) {
    Logger::getInstance()
        .addWarning(
            "PrepareStencilLocalLaxFriedrichs: isDepo does not have enough "
            "elements. Assuming all higher layers are not depo layers.")
        .print();
  }
  // always resize, so it has the correct number of elements
  isDepo.resize(levelSets.size(), false);

  // Begin with the biggest level set (top LS wrapped around all others)
  auto layerIt = levelSets.rbegin();

  // make empty LS which will contain the final top layer
  auto finalTop = SmartPointer<Domain<T, D>>::New(levelSets[0]->getGrid());

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
        auto layerAbove = SmartPointer<Domain<T, D>>::New(*layerIt);
        BooleanOperation<T, D>(layerAbove, *maskIt,
                               BooleanOperationEnum::RELATIVE_COMPLEMENT)
            .apply();

        BooleanOperation<T, D>(finalTop, layerAbove,
                               BooleanOperationEnum::UNION)
            .apply();
        Prune<T, D>(finalTop).apply();
      }
      layerAboveIsDepo = false;
    }
    ++depoIt;
  }

  // If the lowest layer is a depo substrate, add it to the final top layer
  if (layerAboveIsDepo) {
    BooleanOperation<T, D>(finalTop, *layerIt, BooleanOperationEnum::UNION)
        .apply();
  }

  BooleanOperation<T, D>(finalTop, BooleanOperationEnum::INVERT).apply();
  levelSets.back()->deepCopy(finalTop);
}

/// After advection using the SLLF layer wrapping approach is done,
/// restore the original layer wrapping used everywhere else.
template <class T, int D>
void FinalizeStencilLocalLaxFriedrichs(
    std::vector<SmartPointer<Domain<T, D>>> &levelSets) {
  auto layerIt = levelSets.rbegin();
  auto lastIt = ++levelSets.rbegin();

  BooleanOperation<T, D>(*layerIt, BooleanOperationEnum::INVERT).apply();
  BooleanOperation<T, D>(*layerIt, *lastIt, BooleanOperationEnum::UNION)
      .apply();
}

} // namespace viennals
