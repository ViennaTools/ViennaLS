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

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  lsSmartPointer<lsMesh<N>> mesh = nullptr;
  lsSmartPointer<translatorType> translator = nullptr;
  T maxValue = 0.5;
  bool buildTranslator = false;

public:
  lsToDiskMesh() {}

  lsToDiskMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
               lsSmartPointer<lsMesh<N>> passedMesh, T passedMaxValue = 0.5)
      : levelSet(passedLevelSet), mesh(passedMesh), maxValue(passedMaxValue) {}

  lsToDiskMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
               lsSmartPointer<lsMesh<N>> passedMesh,
               lsSmartPointer<translatorType> passedTranslator,
               T passedMaxValue = 0.5)
      : levelSet(passedLevelSet), mesh(passedMesh),
        translator(passedTranslator), maxValue(passedMaxValue) {
    buildTranslator = true;
  }

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setMesh(lsSmartPointer<lsMesh<N>> passedMesh) { mesh = passedMesh; }

  void setTranslator(lsSmartPointer<translatorType> passedTranslator) {
    translator = passedTranslator;
    buildTranslator = true;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsToDiskMesh.")
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

    lsExpand<T, D>(levelSet, (maxValue * 4) + 1).apply();
    lsCalculateNormalVectors<T, D>(levelSet, maxValue).apply();

    const T gridDelta = levelSet->getGrid().getGridDelta();
    const auto &normalVectors =
        *(levelSet->getPointData().getVectorData("Normals"));

    // set up data arrays
    std::vector<N> values;
    std::vector<std::array<N, 3>> normals;

    // save the extent of the resulting mesh
    std::array<N, 3> minimumExtent = {};
    std::array<N, 3> maximumExtent = {};
    for (unsigned i = 0; i < D; ++i) {
      minimumExtent[i] = std::numeric_limits<T>::max();
      maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    values.reserve(normalVectors.size());
    normals.reserve(normalVectors.size());

    const bool buildTranslatorFlag = buildTranslator;
    unsigned long counter = 0;
    if (buildTranslatorFlag) {
      translator->reserve(normalVectors.size());
    }

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || std::abs(it.getValue()) > maxValue) {
        continue;
      }

      unsigned pointId = it.getPointId();

      // insert pointId-counter pair in translator
      if (buildTranslatorFlag) {
        translator->insert({pointId, counter++});
      }

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
        node[i] = double(it.getStartIndices(i)) * gridDelta;

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
      double scaling = it.getValue() * gridDelta * max;
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
      values.push_back(it.getValue());
    }

    mesh->insertNextScalarData(values, "LSValues");
    mesh->insertNextVectorData(normals, "Normals");
    mesh->minimumExtent = minimumExtent;
    mesh->maximumExtent = maximumExtent;
  }
};

#endif // LS_TO_DISK_MESH_HPP
