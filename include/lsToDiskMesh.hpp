#ifndef LS_TO_DISK_MESH_HPP
#define LS_TO_DISK_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>

/// This class creates a mesh from the level set
/// with all grid points with a level set value <= 0.5.
/// These grid points are shifted in space towards the
/// direction of their normal vector by grid delta * LS value.
/// Grid delta and the origin grid point are saved for each point.
/// This allows for a simple setup of disks for ray tracing.
template <class T, int D> class lsToDiskMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsDomain<T, D> *levelSet = nullptr;
  lsMesh *mesh = nullptr;
  T maxValue = 0.5;

public:
  lsToDiskMesh() {}

  lsToDiskMesh(lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh,
               T passedMaxValue = 0.5)
      : levelSet(&passedLevelSet), mesh(&passedMesh), maxValue(passedMaxValue) {
  }

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

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

    mesh->clear();

    lsExpand<T, D>(*levelSet, (maxValue * 4) + 1).apply();
    lsCalculateNormalVectors<T, D>(*levelSet, maxValue).apply();

    const T gridDelta = levelSet->getGrid().getGridDelta();
    const auto &normalVectors = levelSet->getNormalVectors();


    
    // set up data arrays
    //std::vector<double> values(normalVectors.size());
    //std::vector<double> gridSpacing(normalVectors.size());
    //std::vector<std::array<double, 3>> normals(normalVectors.size());
    
    std::vector<double> values;
    std::vector<double> gridSpacing;
    std::vector<std::array<double, 3>> normals;

    values.reserve(normalVectors.size());
    gridSpacing.reserve(normalVectors.size());
    normals.reserve(normalVectors.size());

    int count = 0;

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || std::abs(it.getValue()) > maxValue) {
        continue;
      }

      unsigned pointId = it.getPointId();

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node shifted by ls value in direction of the
      // normal vector
      std::array<double, 3> node;
      node[2] = 0.;
      double max = 0.;
      for (unsigned i = 0; i < D; ++i) {
        // original position
        node[i] = double(it.getStartIndices(i)) * gridDelta;

        if (std::abs(normalVectors[pointId][i]) > max) {
          max = std::abs(normalVectors[pointId][i]);
        }
      }

      // now normalize vector to scale position correctly to manhatten distance
      double scaling = it.getValue() * gridDelta * max;
      for (unsigned i = 0; i < D; ++i) {
        node[i] -= scaling * normalVectors[pointId][i];
      }

      //TODO: REMOVE DEBUG
      if(normalVectors[pointId][0]==0.0 && normalVectors[pointId][1]==0.0 && normalVectors[pointId][2]==0.0)
        std::cout << "blub ";

      mesh->insertNextNode(node);

      // add data into mesh
      //values[pointId] = it.getValue();
      //gridSpacing[pointId] = gridDelta;
      // copy normal
      std::array<double, 3> test_normal;
      if (D == 2)
        normals[pointId][2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        //normals[pointId][i] = normalVectors[pointId][i];
        test_normal[i] = normalVectors[pointId][i];
      }

      normals.push_back(test_normal);
      values.push_back(it.getValue());
      gridSpacing.push_back(gridDelta);


      //TODO: REMOVE
      count++;
    }

    std::cout << "runs trough loop: " << count << std::endl;

    //REVIEW: there are some zero values in the following vectors

    mesh->insertNextScalarData(values, "LSValues");
    mesh->insertNextScalarData(gridSpacing, "gridSpacing");
    mesh->insertNextVectorData(normals, "Normals");
  }
};

#endif // LS_TO_DISK_MESH_HPP
