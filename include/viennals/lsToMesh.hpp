#pragma once

#include <lsPreCompileMacros.hpp>

#include <vector>

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Extract the regular grid, on which the level set values are
/// defined, to an explicit Mesh<>. The Vertices will contain
/// the level set value stored at its location. (This is very useful
/// for debugging)
template <class T, int D> class ToMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  bool onlyDefined = true;
  bool onlyActive = false;
  static constexpr long long maxDomainExtent = 1e6;

public:
  ToMesh() {};

  ToMesh(const SmartPointer<Domain<T, D>> passedLevelSet,
         SmartPointer<Mesh<T>> passedMesh, bool passedOnlyDefined = true,
         bool passedOnlyActive = false)
      : levelSet(passedLevelSet), mesh(passedMesh),
        onlyDefined(passedOnlyDefined), onlyActive(passedOnlyActive) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setOnlyDefined(bool passedOnlyDefined) {
    onlyDefined = passedOnlyDefined;
  }

  void setOnlyActive(bool passedOnlyActive) { onlyActive = passedOnlyActive; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to ToMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance().addWarning("No mesh was passed to ToMesh.").print();
      return;
    }

    mesh->clear();

    // check if level set is empty
    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    // LS data
    std::vector<T> LSValues;
    std::vector<T> subLS;
    // point data
    const auto &pointData = levelSet->getPointData();
    using DomainType = Domain<T, D>;
    using ScalarDataType = typename DomainType::PointDataType::ScalarDataType;
    using VectorDataType = typename DomainType::PointDataType::VectorDataType;

    // get all the data of the LS
    std::vector<ScalarDataType> scalarData;
    std::vector<VectorDataType> vectorData;
    scalarData.resize(pointData.getScalarDataSize());
    vectorData.resize(pointData.getVectorDataSize());

    const T gridDelta = levelSet->getGrid().getGridDelta();

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if ((onlyDefined && !it.isDefined()) ||
          (onlyActive && std::abs(it.getValue()) > 0.5))
        continue;

      if (!onlyDefined && !it.isDefined()) {
        bool skipPoint = false;
        for (unsigned i = 0; i < D; ++i) {
          if (std::abs(it.getStartIndices(i)) > maxDomainExtent) {
            skipPoint = true;
          }
        }
        if (skipPoint) {
          continue;
        }
      }

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node
      std::array<T, 3> node;
      if (D == 2)
        node[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        node[i] = T(it.getStartIndices(i)) * gridDelta;
      }
      mesh->insertNextNode(node);

      // insert LS value
      if (it.isDefined()) {
        LSValues.push_back(it.getDefinedValue());
      } else {
        LSValues.push_back((it.getValue() < 0) ? -1000 : 1000);
      }
      subLS.push_back(it.getSegmentId());

      // add all saved LS data
      for (unsigned i = 0; i < pointData.getScalarDataSize(); ++i) {
        if (const auto dataPointer = pointData.getScalarData(i);
            dataPointer != nullptr) {
          const auto &currentData = *dataPointer;
          scalarData[i].push_back(currentData[it.getPointId()]);
        } else {
          Logger::getInstance()
              .addWarning("ToMesh: Tried to access out of bounds scalarData! "
                          "Ignoring.")
              .print();
          break;
        }
      }

      for (unsigned i = 0; i < pointData.getVectorDataSize(); ++i) {
        if (const auto dataPointer = pointData.getVectorData(i);
            dataPointer != nullptr) {
          const auto &currentData = *dataPointer;
          vectorData[i].push_back(currentData[it.getPointId()]);
        } else {
          Logger::getInstance()
              .addWarning("ToMesh: Tried to access out of bounds vectorData! "
                          "Ignoring.")
              .print();
          break;
        }
      }
    }

    mesh->cellData.insertNextScalarData(LSValues, "LSValues");
    mesh->cellData.insertNextScalarData(subLS, "SegmentID");

    // append all scalar and vector data
    // just move it into the new structure, since we do not need it anymore
    for (unsigned i = 0; i < scalarData.size(); ++i) {
      mesh->cellData.insertNextScalarData(std::move(scalarData[i]),
                                          pointData.getScalarDataLabel(i));
    }

    for (unsigned i = 0; i < vectorData.size(); ++i) {
      mesh->cellData.insertNextVectorData(std::move(vectorData[i]),
                                          pointData.getVectorDataLabel(i));
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToMesh)

} // namespace viennals
