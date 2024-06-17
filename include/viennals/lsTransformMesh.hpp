#pragma once

#include <cmath>

#include <hrleVectorType.hpp>

#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Enumeration for the different types of
/// transformation operations
enum struct TransformEnum : unsigned {
  TRANSLATION = 0,
  ROTATION = 1,
  SCALE = 2
};

template <class T> class TransformMesh {
  SmartPointer<Mesh<T>> mesh = nullptr;
  TransformEnum transform = TransformEnum::TRANSLATION;
  hrleVectorType<double, 3> transformVector{};
  double angle = 0.0;
  double numericEps = 1e-6;

  // check vector for all zeros
  bool isValidVector() {
    if (DotProduct(transformVector, transformVector) < numericEps) {
      return false;
    } else {
      return true;
    }
  }

  void translateMesh() {
    // if vector is 0, just dont do anything
    if (!isValidVector()) {
      return;
    }
    for (auto &node : mesh->nodes) {
      for (unsigned i = 0; i < 3; ++i) {
        node[i] += transformVector[i];
      }
    }
  }

  // rotation of mesh around the vector transform vector by the angle
  void rotateMesh() {
    // if invalid input, dont do anything
    if (!isValidVector()) {
      return;
    }
    if (std::abs(angle) < numericEps) {
      return;
    }
    const double norm = std::sqrt(transformVector[0] * transformVector[0] +
                                  transformVector[1] * transformVector[1] +
                                  transformVector[2] * transformVector[2]);
    for (unsigned i = 0; i < 3; ++i) {
      transformVector[i] /= norm;
    }
    const auto axis = transformVector;
    double sinAngle = std::sin(angle);
    double cosAngle = std::cos(angle);
    double uniMinusCosAngle = 1 - cosAngle;
    for (auto &node : mesh->nodes) {
      auto oldNode = node;
      double term1 = axis[0] * node[0] + axis[1] * node[1] + axis[2] * node[2];
      for (unsigned i = 0; i < 3; ++i) {
        unsigned ip1 = (i + 1) % 3;
        unsigned ip2 = (i + 2) % 3;
        node[i] =
            axis[i] * term1 * uniMinusCosAngle + oldNode[i] * cosAngle +
            (axis[ip1] * oldNode[ip2] - axis[ip2] * oldNode[ip1]) * sinAngle;
      }
    }
  }

  void scaleMesh() {
    if (!isValidVector()) {
      Logger::getInstance()
          .addWarning("TransformMesh: TransformVector is not valid!")
          .print();
      return;
    }
    for (auto &node : mesh->nodes) {
      for (unsigned i = 0; i < 3; ++i) {
        node[i] *= transformVector[i];
      }
    }
  }

public:
  TransformMesh(SmartPointer<Mesh<T>> passedMesh,
                TransformEnum passedTransform = TransformEnum::TRANSLATION,
                std::array<double, 3> passedTransformVector = {},
                double passedAngle = 0.0)
      : mesh(passedMesh), transform(passedTransform),
        transformVector(passedTransformVector), angle(passedAngle) {}

  TransformMesh(SmartPointer<Mesh<T>> passedMesh,
                TransformEnum passedTransform = TransformEnum::TRANSLATION,
                hrleVectorType<double, 3> passedTransformVector = {},
                double passedAngle = 0.0)
      : mesh(passedMesh), transform(passedTransform),
        transformVector(passedTransformVector), angle(passedAngle) {}

  void apply() {
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh passed to TransformMesh. Not transforming!")
          .print();
    }

    switch (transform) {
    case TransformEnum::TRANSLATION:
      translateMesh();
      break;
    case TransformEnum::ROTATION:
      rotateMesh();
      break;
    case TransformEnum::SCALE:
      scaleMesh();
      break;
    default:
      Logger::getInstance()
          .addWarning(
              "Invalid transform passed to TransformMesh. Not transforming!")
          .print();
    }
  }
};

} // namespace viennals
