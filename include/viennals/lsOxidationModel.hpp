#pragma once

#include <lsOxidationMask.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <vcTimer.hpp>

namespace viennals {

template <class T> struct OxidationCouplingParameters {
  unsigned maxIterations = 10;
  T tolerance = 1e-6;
  T relaxation = 1.;
};


/// Iterates diffusion, oxide deformation, and pressure-dependent reaction-rate
/// feedback on the shared Cartesian solve grid.
template <class T, int D> class OxidationModel {
  using IndexType = viennahrle::Index<D>;

  SmartPointer<OxidationDiffusion<T, D>> diffusionField = nullptr;
  SmartPointer<OxidationDeformation<T, D>> deformationField =
      nullptr;
  OxidationCouplingParameters<T> parameters;
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
      SmartPointer<OxidationDeformation<T, D>>
          passedDeformationField,
      OxidationCouplingParameters<T> passedParameters = {})
      : diffusionField(passedDiffusionField),
        deformationField(passedDeformationField), parameters(passedParameters) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationModel>::New(std::forward<Args>(args)...);
  }

  void setDiffusionField(
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField) {
    diffusionField = passedDiffusionField;
  }

  void setDeformationField(
      SmartPointer<OxidationDeformation<T, D>>
          passedDeformationField) {
    deformationField = passedDeformationField;
  }

  void setParameters(OxidationCouplingParameters<T> passedParameters) {
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

    // forEachSolutionNode order is deterministic on a fixed grid, so nodeIndices
    // is collected once and p_raw is refilled in-place on subsequent iterations.
    std::vector<IndexType> nodeIndices;
    std::vector<T>         p_raw;
    std::vector<T>         p_km1, p_k;  // Aitken Δ² history

    for (; iterations < parameters.maxIterations; ++iterations) {
      logDebug("OxidationModel: coupling iteration " +
              std::to_string(iterations + 1) + "/" +
              std::to_string(parameters.maxIterations) +
              " starting diffusion solve");
      diffusionField->apply();
      logDebug("OxidationModel: diffusion solve complete, nodes=" +
              std::to_string(diffusionField->getNumberOfSolutionNodes()) +
              ", iterations=" +
              std::to_string(diffusionField->getIterations()) +
              ", residual=" + std::to_string(diffusionField->getResidual()));
      if (!diffusionField->lastSolveConverged() ||
          !diffusionField->hasFiniteConcentrationField()) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ = "diffusion solve failed (residual=" +
                         std::to_string(diffusionField->getNormalizedResidual()) +
                         ", tolerance=" +
                         std::to_string(diffusionField->getParameters().tolerance) +
                         ")";
        Logger::getInstance()
            .addWarning("OxidationModel: " + failureReason_ + ".")
            .print();
        return;
      }

      logDebug("OxidationModel: coupling iteration " +
              std::to_string(iterations + 1) + "/" +
              std::to_string(parameters.maxIterations) +
              " starting deformation solve");
      Timer<> tDeform;
      tDeform.start();
      deformationField->apply();
      tDeform.finish();
      logDebug("OxidationModel: deformation solve complete, nodes=" +
              std::to_string(deformationField->getNumberOfSolutionNodes()) +
              ", iterations=" +
              std::to_string(deformationField->getIterations()) +
              ", residual=" + std::to_string(deformationField->getResidual()));
      if (!deformationField->lastSolveConverged() ||
          !deformationField->hasFiniteSolution()) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ = "deformation solve failed (mechanics=" +
                         std::to_string(deformationField->getResidual()) +
                         ", pressure=" +
                         std::to_string(deformationField->getLastPressureResidual()) +
                         ", stokes=" +
                         std::to_string(deformationField->getLastStokesResidual()) +
                         ")";
        Logger::getInstance()
            .addWarning("OxidationModel: " + failureReason_ + ".")
            .print();
        return;
      }
      if (Logger::hasTiming())
        Logger::getInstance()
            .addTiming("    deformation(couplingIter=" +
                       std::to_string(iterations + 1) + ")", tDeform)
            .print();

      // Collect raw pressures G(x_k) from deformation.
      if (iterations == 0) {
        nodeIndices.clear(); p_raw.clear();
        deformationField->forEachSolutionNode(
            [&](const IndexType &idx, T p) {
              nodeIndices.push_back(idx);
              p_raw.push_back(p);
            });
      } else {
        std::size_t ii = 0;
        deformationField->forEachSolutionNode(
            [&](const IndexType &, T p) { p_raw[ii++] = p; });
      }
      const std::size_t n = p_raw.size();
      if (!std::all_of(p_raw.begin(), p_raw.end(),
                       [](T value) { return std::isfinite(value); })) {
        residual = std::numeric_limits<T>::infinity();
        failureReason_ = "deformation pressure feedback produced non-finite values";
        Logger::getInstance()
            .addWarning("OxidationModel: " + failureReason_ + ".")
            .print();
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
          const T Fkm1 = p_k[i]   - p_km1[i];
          const T Fk   = p_raw[i] - p_k[i];
          const T dF   = Fk - Fkm1;
          dot_Fkm1_dF += Fkm1 * dF;
          dot_dF_dF   += dF   * dF;
        }
        if (dot_dF_dF > T(1e-100))
          aitkenTheta = std::max(T(0.1),
                        std::min(T(1.0), -dot_Fkm1_dF / dot_dF_dF));
      }

      // Compute blended pressure and relative-change residual.
      std::vector<T> p_blended(n);
      T maxChange = T(0), maxPressure = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        const T oldP  = p_k.empty() ? T(0) : p_k[i];
        p_blended[i]  = (T(1) - aitkenTheta) * oldP + aitkenTheta * p_raw[i];
        maxChange     = std::max(maxChange,   std::abs(p_blended[i] - oldP));
        maxPressure   = std::max(maxPressure, std::abs(p_blended[i]));
      }
      residual = (maxPressure > std::numeric_limits<T>::epsilon())
                 ? maxChange / maxPressure : maxChange;

      // Feed blended pressures to diffusion.
      for (std::size_t i = 0; i < n; ++i)
        diffusionField->setPressure(nodeIndices[i], p_blended[i]);

      // Shift Aitken history; seed p_km1 with zeros on the first iteration.
      p_km1 = p_k.empty() ? std::vector<T>(n, T(0)) : std::move(p_k);
      p_k   = p_blended;

      logDebug("OxidationModel: coupling iteration " +
              std::to_string(iterations + 1) +
              " pressure-feedback residual=" + std::to_string(residual));
      if (residual < parameters.tolerance)
        break;
    }
    const unsigned completedIterations =
        std::min(iterations + 1, parameters.maxIterations);
    logDebug("OxidationModel: coupled solve complete, iterations=" +
            std::to_string(completedIterations) +
            ", residual=" + std::to_string(residual));
    converged_ = std::isfinite(residual) && residual <= parameters.tolerance;
    if (!converged_)
      Logger::getInstance()
          .addWarning("OxidationModel: pressure-concentration coupling did not "
                      "converge after " + std::to_string(completedIterations) +
                      " iterations (residual=" + std::to_string(residual) +
                      ", tolerance=" + std::to_string(parameters.tolerance) +
                      "). Consider increasing maxIterations or relaxation.")
          .print();
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  bool hasConverged() const { return converged_; }
  std::string getFailureReason() const { return failureReason_; }

private:
  static void logDebug(const std::string &message) {
    if (Logger::hasDebug())
      Logger::getInstance().addDebug(message).print();
  }

};

} // namespace viennals
