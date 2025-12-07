#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMaterialMap.hpp>
#include <lsMesh.hpp>
#include <unordered_map>

namespace viennals {

using namespace viennacore;

/// This class creates a mesh from the level set
/// with all grid points with a level set value <= 0.5.
/// These grid points are shifted in space towards the
/// direction of their normal vector by grid delta * LS value.
/// Grid delta and the origin grid point are saved for each point.
/// This allows for a simple setup of disks for ray tracing.
template <class T, int D, class N = T> class ToDiskMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

public:
  using TranslatorType = std::unordered_map<unsigned long, unsigned long>;

private:
  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  SmartPointer<Mesh<N>> mesh = nullptr;
  SmartPointer<TranslatorType> translator = nullptr;
  SmartPointer<MaterialMap> materialMap = nullptr;
  T maxValue = 0.5;
  bool buildTranslator = false;
  static constexpr double wrappingLayerEpsilon = 1e-4;

public:
  ToDiskMesh() = default;

  ToDiskMesh(SmartPointer<Mesh<N>> passedMesh, T passedMaxValue = 0.5)
      : mesh(passedMesh), maxValue(passedMaxValue) {}

  ToDiskMesh(SmartPointer<Domain<T, D>> passedLevelSet,
             SmartPointer<Mesh<N>> passedMesh, T passedMaxValue = 0.5)
      : mesh(passedMesh), maxValue(passedMaxValue) {
    levelSets.push_back(passedLevelSet);
  }

  ToDiskMesh(SmartPointer<Domain<T, D>> passedLevelSet,
             SmartPointer<Mesh<N>> passedMesh,
             SmartPointer<TranslatorType> passedTranslator,
             T passedMaxValue = 0.5)
      : mesh(passedMesh), translator(passedTranslator),
        maxValue(passedMaxValue) {
    levelSets.push_back(passedLevelSet);
    buildTranslator = true;
  }

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  /// Pushes the passed level set to the back of the list of level sets
  void insertNextLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void setMesh(SmartPointer<Mesh<N>> passedMesh) { mesh = passedMesh; }

  void setTranslator(SmartPointer<TranslatorType> passedTranslator) {
    translator = passedTranslator;
    buildTranslator = true;
  }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void clearLevelSets() { levelSets.clear(); }

  void apply() {
    if (levelSets.empty()) {
      Logger::getInstance()
          .addError("No level sets passed to ToDiskMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToDiskMesh.")
          .print();
      return;
    }
    if (buildTranslator && translator == nullptr) {
      VIENNACORE_LOG_WARNING("No translator was passed to ToDiskMesh.");
    }

    mesh->clear();

    // expand top levelset
    Expand<T, D>(levelSets.back(), (maxValue * 4) + 1).apply();
    CalculateNormalVectors<T, D>(levelSets.back(), maxValue).apply();

    const T gridDelta = levelSets.back()->getGrid().getGridDelta();
    const auto &normalVectors =
        *(levelSets.back()->getPointData().getVectorData(
            CalculateNormalVectors<T, D>::normalVectorsLabel));

    // set up data arrays
    std::vector<N> values;
    std::vector<Vec3D<N>> normals;
    std::vector<N> materialIds;

    // save the extent of the resulting mesh
    Vec3D<N> minimumExtent;
    Vec3D<N> maximumExtent;
    for (unsigned i = 0; i < D; ++i) {
      minimumExtent[i] = std::numeric_limits<T>::max();
      maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    values.reserve(normalVectors.size());
    normals.reserve(normalVectors.size());
    materialIds.reserve(normalVectors.size());

    const bool useMaterialMap = materialMap != nullptr;
    const bool buildTranslatorFlag = buildTranslator;
    unsigned long counter = 0;
    if (buildTranslatorFlag) {
      translator->clear();
      translator->reserve(normalVectors.size());
    }

    // an iterator for each levelset
    std::vector<viennahrle::ConstSparseIterator<hrleDomainType>> iterators;
    for (const auto levelSet : levelSets) {
      iterators.emplace_back(levelSet->getDomain());
    }

    // iterate over top levelset
    for (auto &topIt = iterators.back(); !topIt.isFinished(); ++topIt) {
      if (!topIt.isDefined() || std::abs(topIt.getValue()) > maxValue) {
        continue;
      }

      unsigned pointId = topIt.getPointId();

      // insert pointId-counter pair in translator
      if (buildTranslatorFlag) {
        translator->insert({pointId, counter++});
      }

      // insert material ID
      const T value = topIt.getValue();
      int matId = levelSets.size() - 1;
      for (int lowerLevelSetId = 0; lowerLevelSetId < levelSets.size() - 1;
           ++lowerLevelSetId) {
        // check if there is any other levelset at the same point:
        // put iterator to same position as the top levelset
        iterators[lowerLevelSetId].goToIndicesSequential(
            topIt.getStartIndices());
        if (iterators[lowerLevelSetId].getValue() <=
            value + wrappingLayerEpsilon) {
          matId = lowerLevelSetId;
          break;
        }
      }
      if (useMaterialMap)
        matId = materialMap->getMaterialId(matId);

      materialIds.push_back(matId);

      // insert vertex
      std::array<unsigned, 1> vertex{};
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node shifted by ls value in direction of the
      // normal vector
      Vec3D<N> node;
      node[2] = 0.;
      double max = 0.;
      for (unsigned i = 0; i < D; ++i) {
        // original position
        node[i] = double(topIt.getStartIndices(i)) * gridDelta;

        // save extent
        if (node[i] < minimumExtent[i]) {
          minimumExtent[i] = node[i];
        } else if (node[i] > maximumExtent[i]) {
          maximumExtent[i] = node[i];
        }

        if (std::abs(normalVectors[pointId][i]) > max) {
          max = std::abs(normalVectors[pointId][i]);
        }
      }

      // now normalize vector to scale position correctly to manhatten distance
      double scaling = value * gridDelta * max;
      for (unsigned i = 0; i < D; ++i) {
        node[i] -= scaling * normalVectors[pointId][i];
      }

      mesh->insertNextNode(node);

      // add data into mesh
      // copy normal
      Vec3D<N> normal;
      if (D == 2)
        normal[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        normal[i] = normalVectors[pointId][i];
      }

      normals.push_back(normal);
      values.push_back(value);
    }

    // delete normal vectors point data again as it is not needed anymore
    {
      auto &pointData = levelSets.back()->getPointData();
      auto index = pointData.getVectorDataIndex(
          CalculateNormalVectors<T, D>::normalVectorsLabel);
      if (index < 0) {
        VIENNACORE_LOG_WARNING(
            "ToDiskMesh: Could not find normal vector data.");
      } else {
        pointData.eraseVectorData(index);
      }
    }

    mesh->cellData.insertNextScalarData(values, "LSValues");
    mesh->cellData.insertNextVectorData(normals, "Normals");
    mesh->cellData.insertNextScalarData(materialIds, "MaterialIds");
    mesh->minimumExtent = minimumExtent;
    mesh->maximumExtent = maximumExtent;
  }
};

PRECOMPILE_PRECISION_DIMENSION(ToDiskMesh)

} // namespace viennals
