#pragma once

#include <hrleSparseBoxIterator.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>
#include <lsVelocityField.hpp>

#include <vcVectorType.hpp>

// Choose adaptive epsilon strategy
// 0: Original fixed epsilon
// 1: Simple adaptive epsilon (recommended)
// 2: Full adaptive epsilon (most accurate but slower)
#ifndef ADAPTIVE_EPSILON_STRATEGY
#define ADAPTIVE_EPSILON_STRATEGY 0
#endif

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
  const double baseNormalEpsilon =
      std::cbrt(std::numeric_limits<double>::epsilon());

  // Final dissipation coefficients that are used by the time integrator. If
  // D==2 last entries are 0.
  Vec3D<T> finalAlphas;
  static constexpr unsigned numStencilPoints = hrleUtil::pow(2 * order + 1, D);
  static double maxDissipation; // default: std::numeric_limits<double>::max();

  // Adaptive epsilon control parameters
  static constexpr T velocityScaleFactor =
      0.5; // Controls velocity-based scaling
  static constexpr T smoothnessThreshold =
      0.1; // Threshold for smoothness detection
  static constexpr T epsilonMinScale = 0.1;   // Minimum epsilon scale factor
  static constexpr T epsilonMaxScale = 100.0; // Maximum epsilon scale factor

  static T pow2(const T &value) { return value * value; }

  // Calculate adaptive epsilon based on velocity magnitude, grid resolution,
  // and local velocity smoothness
  T calculateAdaptiveEpsilon(T velocityMagnitude, T gridDelta,
                             const Vec3D<T> &localCoordArray, int material,
                             const Vec3D<T> &normal) const {
    // Base epsilon scaled by grid resolution
    T adaptiveEpsilon = baseNormalEpsilon * gridDelta;

    // Scale with velocity magnitude to maintain numerical stability
    // Use larger epsilon for larger velocities to avoid catastrophic
    // cancellation
    T absVelocity = std::abs(velocityMagnitude);
    if (absVelocity > 1e-12) {
      T velocityScale = 1.0 + velocityScaleFactor * std::log(1.0 + absVelocity);
      adaptiveEpsilon *= velocityScale;
    }

    // Estimate local velocity smoothness using a simplified approach
    // Only check smoothness if velocity is significant
    if (absVelocity > 1e-10) {
      T maxRelativeVariation = 0.0;

      // Sample velocity in each coordinate direction with small perturbation
      for (int dir = 0; dir < D; ++dir) {
        Vec3D<T> perturbedNormal = normal;
        perturbedNormal[dir] += smoothnessThreshold;

        // Normalize the perturbed normal
        T normSq = 0.0;
        for (int i = 0; i < D; ++i) {
          normSq += perturbedNormal[i] * perturbedNormal[i];
        }
        if (normSq > 1e-12) {
          T invNorm = 1.0 / std::sqrt(normSq);
          for (int i = 0; i < D; ++i) {
            perturbedNormal[i] *= invNorm;
          }

          T perturbedVel = velocities->getScalarVelocity(
              localCoordArray, material, perturbedNormal,
              neighborIterator.getCenter().getPointId());

          T relativeVariation =
              std::abs(perturbedVel - velocityMagnitude) / absVelocity;
          maxRelativeVariation =
              std::max(maxRelativeVariation, relativeVariation);
        }
      }

      // If velocity varies significantly, use larger epsilon for stability
      if (maxRelativeVariation > smoothnessThreshold) {
        T smoothnessFactor = 1.0 + maxRelativeVariation;
        adaptiveEpsilon *= smoothnessFactor;
      }
    }

    // Clamp epsilon to reasonable bounds
    T minEpsilon = baseNormalEpsilon * gridDelta * epsilonMinScale;
    T maxEpsilon = baseNormalEpsilon * gridDelta * epsilonMaxScale;

    return std::max(minEpsilon, std::min(maxEpsilon, adaptiveEpsilon));
  }

  Vec3D<T> calculateNormal(const viennahrle::Index<D> &offset) {
    Vec3D<T> normal = {0.0, 0.0, 0.0};
    constexpr int startIndex = -1;
    T modulus = 0.;

    for (unsigned i = 0; i < D; ++i) {
      viennahrle::Index<D> index(offset);
      std::array<T, 3> values;
      for (unsigned j = 0; j < 3; ++j) {
        index[i] = startIndex + j;
        values[j] = neighborIterator.getNeighbor(index).getValue();
      }
      normal[i] = FiniteDifferences<T>::calculateGradient(
          values.data(), levelSet->getGrid().getGridDelta());
      modulus += normal[i] * normal[i];
    }
    modulus = 1 / std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      normal[i] *= modulus;
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
      std::array<T, numValues> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values[j] = neighborIterator.getNeighbor(index).getValue();
      }

      gradient[i] =
          FiniteDifferences<T, finiteDifferenceScheme>::calculateGradient(
              values.data(), levelSet->getGrid().getGridDelta());
    }

    return gradient;
  }

  VectorType<T, D> calculateGradientDiff() {
    VectorType<T, D> gradient;

    constexpr unsigned numValues =
        FiniteDifferences<T, finiteDifferenceScheme>::getNumberOfValues();
    const int startIndex =
        -std::floor(numValues / 2); // std::floor constexpr in C++23

    for (unsigned i = 0; i < D; ++i) {
      viennahrle::Index<D> index(0);
      std::array<T, numValues> values;
      for (unsigned j = 0; j < numValues; ++j) {
        index[i] = startIndex + j;
        values[j] = neighborIterator.getNeighbor(index).getValue();
      }

      gradient[i] =
          FiniteDifferences<T, finiteDifferenceScheme>::calculateGradientDiff(
              values.data(), levelSet->getGrid().getGridDelta());
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

  // Simple adaptive epsilon calculation for better performance
  T calculateSimpleAdaptiveEpsilon(T velocityMagnitude, T gridDelta) const {
    T absVel = std::abs(velocityMagnitude);
    T baseEps = baseNormalEpsilon * gridDelta;

    // Simple velocity-based scaling to avoid numerical issues
    if (absVel > 1e-10) {
      // Scale epsilon with velocity magnitude (clamped to reasonable range)
      T scale = std::max(T(0.1), std::min(T(10.0), absVel));
      return baseEps * scale;
    }

    return baseEps;
  }

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
    denominator = 1 / std::sqrt(denominator);
    for (unsigned i = 0; i < D; ++i) {
      normalVector[i] *= denominator;
    }

    double scalarVelocity = velocities->getScalarVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());
    auto vectorVelocity = velocities->getVectorVelocity(
        coordinate, material, normalVector,
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
      std::array<VectorType<T, D>, numStencilPoints> alphas{};

      viennahrle::Index<D> currentIndex(-order);
      for (size_t i = 0; i < numStencilPoints; ++i) {
        VectorType<T, D> alpha;
        Vec3D<T> normal = calculateNormal(currentIndex);

        // Check for corrupted normal
        if ((std::abs(normal[0]) < 1e-6) && (std::abs(normal[1]) < 1e-6) &&
            (std::abs(normal[2]) < 1e-6)) {
          alphas[i].fill(0);
          continue;
        }

        Vec3D<T> normal_p{normal[0], normal[1], normal[2]};
        Vec3D<T> normal_n{normal[0], normal[1], normal[2]};

        VectorType<T, D> velocityDelta;
        velocityDelta.fill(0.);

        // get local velocity
        Vec3D<T> localCoordArray = coordinate;
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

        // Calculate epsilon based on selected strategy
        T DN;
#if ADAPTIVE_EPSILON_STRATEGY == 0
        // Original fixed epsilon approach
        DN = std::abs(baseNormalEpsilon * scalarVelocity);
#elif ADAPTIVE_EPSILON_STRATEGY == 1
        // Simple adaptive epsilon (recommended)
        DN = calculateSimpleAdaptiveEpsilon(localScalarVelocity,
                                            levelSet->getGrid().getGridDelta());
#elif ADAPTIVE_EPSILON_STRATEGY == 2
        // Full adaptive epsilon (most accurate but slower)
        DN = calculateAdaptiveEpsilon(localScalarVelocity,
                                      levelSet->getGrid().getGridDelta(),
                                      localCoordArray, material, normal);
#else
        // Fallback to original method
        DN = std::abs(baseNormalEpsilon * scalarVelocity);
#endif

        for (int k = 0; k < D; ++k) {

          normal_p[k] -= DN; // p=previous
          normal_n[k] += DN; // n=next

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
            alpha[k] = maxDissipation;
          }
        }

        alphas[i] = alpha;

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
