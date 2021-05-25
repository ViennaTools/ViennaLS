#ifndef LS_TO_DISK_MESH_HPP
#define LS_TO_DISK_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>
#include <unordered_map>

/// This class creates a mesh from the level set
/// with all grid points with a level set value <= 0.5.
/// These grid points are shifted in space towards the
/// direction of their normal vector by grid delta * LS value.
/// Grid delta and the origin grid point are saved for each point.
/// This allows for a simple setup of disks for ray tracing.
template <class T, int D, class N = T> class lsToDiskMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;
  typedef std::unordered_map<unsigned long, unsigned long> translatorType;

  std::vector<lsSmartPointer<lsDomain<T, D>>> levelSets;
  lsSmartPointer<lsMesh<N>> mesh = nullptr;
  lsSmartPointer<translatorType> translator = nullptr;
  T maxValue = 0.5;
  bool buildTranslator = false;
  static constexpr double wrappingLayerEpsilon = 1e-4;

public:
  lsToDiskMesh() {}

  lsToDiskMesh(lsSmartPointer<lsMesh<N>> passedMesh, T passedMaxValue = 0.5)
      : mesh(passedMesh), maxValue(passedMaxValue) {}

  lsToDiskMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
               lsSmartPointer<lsMesh<N>> passedMesh, T passedMaxValue = 0.5)
      : mesh(passedMesh), maxValue(passedMaxValue) {
    levelSets.push_back(passedLevelSet);
  }

  lsToDiskMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
               lsSmartPointer<lsMesh<N>> passedMesh,
               lsSmartPointer<translatorType> passedTranslator,
               T passedMaxValue = 0.5)
      : mesh(passedMesh), translator(passedTranslator),
        maxValue(passedMaxValue) {
    levelSets.push_back(passedLevelSet);
    buildTranslator = true;
  }

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  /// Pushes the passed level set to the back of the list of level sets
  void insertNextLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void setMesh(lsSmartPointer<lsMesh<N>> passedMesh) { mesh = passedMesh; }

  void setTranslator(lsSmartPointer<translatorType> passedTranslator) {
    translator = passedTranslator;
    buildTranslator = true;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSets.size() < 1) {
      lsMessage::getInstance()
          .addWarning("No level sets passed to lsToDiskMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToDiskMesh.")
          .print();
      return;
    }
    if (buildTranslator && translator == nullptr) {
      lsMessage::getInstance()
          .addWarning("No translator was passed to lsToDiskMesh.")
          .print();
    }

    mesh->clear();

    // expand top levelset
    lsExpand<T, D>(levelSets.back(), (maxValue * 4) + 1).apply();
    lsCalculateNormalVectors<T, D>(levelSets.back(), maxValue).apply();

    const T gridDelta = levelSets.back()->getGrid().getGridDelta();
    const auto &normalVectors =
        *(levelSets.back()->getPointData().getVectorData("Normals"));

    // set up data arrays
    std::vector<N> values;
    std::vector<std::array<N, 3>> normals;
    std::vector<N> materialIds;

    // save the extent of the resulting mesh
    std::array<N, 3> minimumExtent = {};
    std::array<N, 3> maximumExtent = {};
    for (unsigned i = 0; i < D; ++i) {
      minimumExtent[i] = std::numeric_limits<T>::max();
      maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    values.reserve(normalVectors.size());
    normals.reserve(normalVectors.size());
    materialIds.reserve(normalVectors.size());

    const bool buildTranslatorFlag = buildTranslator;
    unsigned long counter = 0;
    if (buildTranslatorFlag) {
      translator->clear();
      translator->reserve(normalVectors.size());
    }

    // an iterator for each levelset
    std::vector<hrleConstSparseIterator<hrleDomainType>> iterators;
    for (const auto levelSet : levelSets) {
      iterators.push_back(
          hrleConstSparseIterator<hrleDomainType>(levelSet->getDomain()));
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
      materialIds.push_back(matId);

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node shifted by ls value in direction of the
      // normal vector
      std::array<N, 3> node;
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
      std::array<N, 3> normal;
      if (D == 2)
        normal[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        normal[i] = normalVectors[pointId][i];
      }

      normals.push_back(normal);
      values.push_back(value);
    }

    mesh->insertNextScalarData(values, "LSValues");
    mesh->insertNextVectorData(normals, "Normals");
    mesh->insertNextScalarData(materialIds, "MaterialIds");
    mesh->minimumExtent = minimumExtent;
    mesh->maximumExtent = maximumExtent;
  }
};

#endif // LS_TO_DISK_MESH_HPP
