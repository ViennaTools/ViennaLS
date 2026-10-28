#pragma once

#include <lsOxidationMask.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include <vcTimer.hpp>

namespace viennals {

struct OxidationCouplingParameters {
  unsigned maxIterations = 10;
  double tolerance = 1e-6;
  double relaxation = 1.;
};

/// Iterates diffusion, oxide deformation, and the stress feedback (normal
/// stress on the reaction rate, pressure on the diffusivity) on the shared
/// Cartesian solve grid.
template <class T, int D> class OxidationModel {
  using IndexType = viennahrle::Index<D>;

  SmartPointer<OxidationDiffusion<T, D>> diffusionField = nullptr;
  SmartPointer<OxidationDeformation<T, D>> deformationField = nullptr;
  OxidationCouplingParameters parameters;
  IndexType minIndex{};
  IndexType maxIndex{};
  bool useRequestedBounds = false;
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  bool converged_ = false;
  std::string failureReason_;

public:
  OxidationModel() = default;

  OxidationModel(
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField,
      SmartPointer<OxidationDeformation<T, D>> passedDeformationField,
      OxidationCouplingParameters passedParameters = {})
      : diffusionField(passedDiffusionField),
        deformationField(passedDeformationField), parameters(passedParameters) {
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationModel>::New(std::forward<Args>(args)...);
  }

  void setDiffusionField(
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField) {
    diffusionField = passedDiffusionField;
  }

  void setDeformationField(
      SmartPointer<OxidationDeformation<T, D>> passedDeformationField) {
    deformationField = passedDeformationField;
  }

  void setParameters(OxidationCouplingParameters passedParameters) {
    parameters = passedParameters;
  }

  void setSolveBounds(const IndexType &passedMinIndex,
                      const IndexType &passedMaxIndex) {
    minIndex = passedMinIndex;
    maxIndex = passedMaxIndex;
    useRequestedBounds = true;
  }

  void clearSolveBounds() { useRequestedBounds = false; }

  void apply() {
    iterations = 0;
    residual = 0.;
    converged_ = false;
    failureReason_.clear();
    if (diffusionField == nullptr || deformationField == nullptr) {
      Logger::getInstance()
          .addError("OxidationModel: Missing diffusion or deformation "
                    "field.")
          .print();
      return;
    }

    if (useRequestedBounds) {
      diffusionField->setSolveBounds(minIndex, maxIndex);
      deformationField->setSolveBounds(minIndex, maxIndex);
    } else {
      diffusionField->clearSolveBounds();
      deformationField->clearSolveBounds();
    }

    // forEachSolutionNode order is deterministic on a fixed grid, so
    // nodeIndices is collected once and the feedback vector is refilled
    // in-place on subsequent iterations.  Per node it holds the pressure
    // and, with a reaction activation volume, the 9 Cauchy stress
    // components (stride 10).
    const auto oxParams = diffusionField->getParameters();
    const bool stressFeedback = oxParams.reactionActivationVolume != 0.;
    const bool pressureFeedback = oxParams.diffusionActivationVolume != 0.;
    const T thermal = T(1.380649e-23) * std::max(T(oxParams.temperature), T(1));
    const std::size_t stride = stressFeedback ? 10 : 1;
    std::vector<IndexType> nodeIndices;
    std::vector<T> p_raw;
    std::vector<T> p_km1, p_k; // Aitken Δ² history
    // One viscosity pass per coupling iteration (setCouplingIterate), off
    // again on every exit of this solve.
    struct CouplingIterateReset {
      SmartPointer<OxidationDeformation<T, D>> field;
      ~CouplingIterateReset() { field->setCouplingIterate(false); }
    } couplingIterateReset{deformationField};

    for (; iterations < parameters.maxIterations; ++iterations) {
      VIENNACORE_LOG_DEBUG("OxidationModel: coupling iteration " +
                           std::to_string(iterations + 1) + "/" +
                           std::to_string(parameters.maxIterations) +
                           " starting diffusion solve");
      diffusionField->apply();
      VIENNACORE_LOG_DEBUG(
          "OxidationModel: diffusion solve complete, nodes=" +
          std::to_string(diffusionField->getNumberOfSolutionNodes()) +
          ", iterations=" + std::to_string(diffusionField->getIterations()) +
          ", residual=" + std::to_string(diffusionField->getResidual()));
      if (!diffusionField->lastSolveConverged() ||
          !diffusionField->hasFiniteConcentrationField()) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ =
            "diffusion solve failed (residual=" +
            std::to_string(diffusionField->getNormalizedResidual()) +
            ", tolerance=" +
            std::to_string(diffusionField->getParameters().tolerance) + ")";
        Logger::getInstance()
            .addWarning("OxidationModel: " + failureReason_ + ".")
            .print();
        return;
      }

      VIENNACORE_LOG_DEBUG("OxidationModel: coupling iteration " +
                           std::to_string(iterations + 1) + "/" +
                           std::to_string(parameters.maxIterations) +
                           " starting deformation solve");
      if (stressFeedback || pressureFeedback)
        deformationField->setCouplingIterate(true);
      Timer<> tDeform;
      tDeform.start();
      deformationField->apply();
      tDeform.finish();
      VIENNACORE_LOG_DEBUG(
          "OxidationModel: deformation solve complete, nodes=" +
          std::to_string(deformationField->getNumberOfSolutionNodes()) +
          ", iterations=" + std::to_string(deformationField->getIterations()) +
          ", residual=" + std::to_string(deformationField->getResidual()));
      if (!deformationField->lastSolveConverged() ||
          !deformationField->hasFiniteSolution()) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ =
            "deformation solve failed (mechanics=" +
            std::to_string(deformationField->getResidual()) + ", continuity=" +
            std::to_string(deformationField->getLastContinuityResidual()) +
            ", stokes=" +
            std::to_string(deformationField->getLastStokesResidual()) + ")";
        VIENNACORE_LOG_WARNING("OxidationModel: " + failureReason_ + ".");
        return;
      }
      if (Logger::hasTiming())
        Logger::getInstance()
            .addTiming("    deformation(couplingIter=" +
                           std::to_string(iterations + 1) + ")",
                       tDeform)
            .print();

      // Collect the raw feedback G(x_k) from deformation.
      if (iterations == 0) {
        nodeIndices.clear();
        p_raw.clear();
        deformationField->forEachSolutionNode(
            [&](const IndexType &idx, T p, const std::array<T, 9> &stress) {
              nodeIndices.push_back(idx);
              p_raw.push_back(p);
              if (stressFeedback)
                p_raw.insert(p_raw.end(), stress.begin(), stress.end());
            });
      } else {
        std::size_t ii = 0;
        deformationField->forEachSolutionNode(
            [&](const IndexType &, T p, const std::array<T, 9> &stress) {
              p_raw[ii++] = p;
              if (stressFeedback)
                for (const T component : stress)
                  p_raw[ii++] = component;
            });
      }
      const std::size_t n = p_raw.size();
      // The first iterate is measured against the feedback the diffusion
      // solve actually used (carried over from the previous step), not
      // against zero: measured against zero, any stress below tolerance*kT/V
      // passed after one iteration, so the stress->rate update ran once per
      // time step — an explicit lag that let a stiff ripple double every step.
      if (iterations == 0 && p_k.empty()) {
        p_k.reserve(n);
        for (const auto &idx : nodeIndices) {
          T pressure;
          std::array<T, 9> stress;
          diffusionField->getFeedback(idx, pressure, stress);
          p_k.push_back(pressure);
          if (stressFeedback)
            p_k.insert(p_k.end(), stress.begin(), stress.end());
        }
      }
      if (!std::all_of(p_raw.begin(), p_raw.end(),
                       [](T value) { return std::isfinite(value); })) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ =
            "deformation pressure feedback produced non-finite values";
        VIENNACORE_LOG_WARNING("OxidationModel: " + failureReason_ + ".");
        return;
      }

      // Aitken Δ²: find θ minimising ||(1-θ)F_{k-1} + θF_k||².
      // Clamped to [0.1, 1.0]: no extrapolation (θ>1) because the
      // pressure→concentration feedback is nonlinear (exponential kinetics)
      // and overshooting the pressure causes oscillatory divergence.
      T aitkenTheta = T(1);
      if (p_k.size() == n && p_km1.size() == n) {
        T dot_Fkm1_dF = T(0), dot_dF_dF = T(0);
        for (std::size_t i = 0; i < n; ++i) {
          const T Fkm1 = p_k[i] - p_km1[i];
          const T Fk = p_raw[i] - p_k[i];
          const T dF = Fk - Fkm1;
          dot_Fkm1_dF += Fkm1 * dF;
          dot_dF_dF += dF * dF;
        }
        if (dot_dF_dF > T(1e-100))
          aitkenTheta =
              std::max(T(0.1), std::min(T(1.0), -dot_Fkm1_dF / dot_dF_dF));
      }

      // Blend, and measure convergence on what the laws actually use: the
      // pressure at every node (diffusivity law, V_D) and the interface
      // normal stress after the median + smoothing (reaction law, V_k).  Raw
      // tensor components are the wrong measure: nodes the moving Si front
      // has just uncovered carry no previous stress, so a max over them read
      // 1.0 every step although the law fills them from their neighbours.
      // Each change is relative to its own magnitude, floored at kT/V (a
      // change small against kT/V does not move the rate).
      std::vector<T> p_blended(n);
      T maxChange = T(0), maxPressure = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        const T oldP = p_k.empty() ? T(0) : p_k[i];
        p_blended[i] = (T(1) - aitkenTheta) * oldP + aitkenTheta * p_raw[i];
        if (i % stride == 0) {
          maxChange = std::max(maxChange, std::abs(p_blended[i] - oldP));
          maxPressure = std::max(maxPressure, std::abs(p_blended[i]));
        }
      }
      residual = T(0);
      if (pressureFeedback)
        residual =
            maxChange /
            std::max(maxPressure,
                     thermal / std::abs(T(oxParams.diffusionActivationVolume)));
      if (stressFeedback) {
        typename OxidationDiffusion<T, D>::StressLookup candidate;
        candidate.reserve(nodeIndices.size());
        for (std::size_t i = 0; i < nodeIndices.size(); ++i) {
          std::array<T, 9> stress;
          std::copy_n(p_blended.begin() + i * stride + 1, 9, stress.begin());
          candidate[detail::gridIndexHash<D>(nodeIndices[i])] = stress;
        }
        const auto next = diffusionField->interfaceNormalStressFor(candidate);
        const auto &used = diffusionField->getInterfaceNormalStress();
        const T unstressed = -T(oxParams.referencePressure);
        T maxStressChange = T(0), maxStress = T(0);
        for (std::size_t i = 0; i < next.size() && i < used.size(); ++i) {
          if (!std::isfinite(next[i]) && !std::isfinite(used[i]))
            continue;
          const T a = std::isfinite(next[i]) ? next[i] : unstressed;
          const T b = std::isfinite(used[i]) ? used[i] : unstressed;
          maxStressChange = std::max(maxStressChange, std::abs(a - b));
          maxStress = std::max(maxStress, std::abs(a));
        }
        residual = std::max(
            residual,
            maxStressChange /
                std::max(maxStress,
                         thermal /
                             std::abs(T(oxParams.reactionActivationVolume))));
      }

      // Converged: keep the feedback the last solves used.  lsOxidation marks
      // the diffusion solved before advecting (no re-solve), so feeding the
      // new values here made the Si front move with a rate that neither the
      // concentration nor the oxide flow had been solved with — an explicit
      // lag that let a stiff stress ripple double every step.
      if (residual < parameters.tolerance) {
        // The advection uses the velocities of this last solve: finish its
        // viscosity iteration (warm-started: a few iterations).
        deformationField->setCouplingIterate(false);
        if (stressFeedback || pressureFeedback)
          deformationField->apply();
        break;
      }

      // Feed the blended pressures and stresses to diffusion.
      for (std::size_t i = 0; i < nodeIndices.size(); ++i) {
        diffusionField->setPressure(nodeIndices[i], p_blended[i * stride]);
        if (stressFeedback) {
          std::array<T, 9> stress;
          std::copy_n(p_blended.begin() + i * stride + 1, 9, stress.begin());
          diffusionField->setStress(nodeIndices[i], stress);
        }
      }

      // Shift Aitken history; seed p_km1 with zeros on the first iteration.
      p_km1 = p_k.empty() ? std::vector<T>(n, T(0)) : std::move(p_k);
      p_k = p_blended;

      VIENNACORE_LOG_DEBUG(
          "OxidationModel: coupling iteration " +
          std::to_string(iterations + 1) +
          " pressure-feedback residual=" + std::to_string(residual));
      if (residual < parameters.tolerance)
        break;
    }
    const unsigned completedIterations =
        std::min(iterations + 1, parameters.maxIterations);
    VIENNACORE_LOG_DEBUG("OxidationModel: coupled solve complete, iterations=" +
                         std::to_string(completedIterations) +
                         ", residual=" + std::to_string(residual));
    converged_ = std::isfinite(residual) && residual <= parameters.tolerance;
    if (!converged_)
      VIENNACORE_LOG_WARNING(
          "OxidationModel: pressure-concentration coupling did not "
          "converge after " +
          std::to_string(completedIterations) +
          " iterations (residual=" + std::to_string(residual) +
          ", tolerance=" + std::to_string(parameters.tolerance) +
          "). Consider increasing maxIterations or relaxation.");
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  bool hasConverged() const { return converged_; }
  std::string getFailureReason() const { return failureReason_; }
};

} // namespace viennals
