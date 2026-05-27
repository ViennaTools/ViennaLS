#pragma once

#include <lsOxidationMask.hpp>

#include <unordered_map>

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

    std::unordered_map<std::size_t, T> previousPressures;
    for (; iterations < parameters.maxIterations; ++iterations) {
      diffusionField->apply();
      deformationField->apply();

      residual = updateDiffusionPressureField(previousPressures);
      if (residual < parameters.tolerance)
        break;
    }

    diffusionField->apply();
    deformationField->apply();
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }

private:
  T updateDiffusionPressureField(
      std::unordered_map<std::size_t, T> &previousPressures) {
    T maxChange = 0.;
    T maxPressure = 0.;
    deformationField->forEachSolutionNode(
        [&](const IndexType &index, T newPressure) {
      const auto key = detail::gridIndexHash<D>(index);
      const auto previousIt = previousPressures.find(key);
      const T oldPressure =
          (previousIt == previousPressures.end()) ? T(0) : previousIt->second;
      const T relaxedPressure =
          parameters.relaxation * newPressure +
          (T(1) - parameters.relaxation) * oldPressure;

      previousPressures[key] = relaxedPressure;
      diffusionField->setPressure(index, relaxedPressure);
      maxChange = std::max(maxChange, std::abs(relaxedPressure - oldPressure));
      maxPressure = std::max(maxPressure, std::abs(relaxedPressure));
        });
    if (maxPressure <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxPressure;
  }
};

} // namespace viennals
