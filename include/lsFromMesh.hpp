#ifndef LS_FROM_MESH_HPP
#define LS_FROM_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <lsDomain.hpp>
#include <lsMesh.hpp>

/// Import the regular grid, on which the level set values are
/// defined, from an explicit lsMesh<>. The Vertices must be defined,
/// as well as a scalar data field "LSValues". If used for custom
/// read-in, make sure all vertices are lexicographically sorted.
template <class T, int D> class lsFromMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  bool sortPointList = true;

public:
  lsFromMesh(){};

  lsFromMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
             const lsSmartPointer<lsMesh<T>> passedMesh)
      : levelSet(passedLevelSet), mesh(passedMesh) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(const lsSmartPointer<lsMesh<T>> passedMesh) {
    mesh = passedMesh;
  }

  void setSortPointList(bool passedSortPointList) {
    sortPointList = passedSortPointList;
  }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsFromMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No lsMesh<> was supplied to lsFromMesh.")
          .print();
      return;
    }

    auto &domain = levelSet->getDomain();
    auto &nodes = mesh->getNodes();
    auto values = mesh->getScalarData("LSValues");

    if (values == nullptr) {
      lsMessage::getInstance()
          .addWarning("Mesh does not contain level set values (\"LSValues\").")
          .print();
      return;
    }

    domain.initialize();

    // if there are no points, just initialize an empty hrleDomain
    if (nodes.empty()) {
      return;
    }

    const hrleGrid<D> &grid = domain.getGrid();
    const T gridDelta = grid.getGridDelta();

    if (hrleVectorType<T, D>(nodes.front()) != grid.getMinGridPoint()) {
      domain.insertNextUndefinedPoint(0, grid.getMinGridPoint(),
                                      (values->front() < 0)
                                          ? lsDomain<T, D>::NEG_VALUE
                                          : lsDomain<T, D>::POS_VALUE);
    }

    hrleVectorType<hrleIndexType, D> lastIndex(nodes.front());
    for (unsigned i = 0; i < D; ++i) {
      lastIndex[i] = nodes.front()[i] / gridDelta;
    }

    auto pointDataIt = nodes.begin();
    auto pointDataEnd = nodes.end();

    auto valueIt = values->begin();
    auto valueEnd = values->end();

    hrleVectorType<bool, D> signs(values->front() < 0);

    while (pointDataIt != pointDataEnd && valueIt != valueEnd) {
      // only read in points within the first 5 layers, to ignore
      // undefined points
      if (std::abs(*valueIt) > 2.5) {
        ++pointDataIt;
        ++valueIt;
        continue;
      }

      bool setPoint = true;

      hrleVectorType<hrleIndexType, D> currentIndex;
      for (unsigned i = 0; i < D; ++i) {
        currentIndex[i] = std::round(pointDataIt->at(i) / gridDelta);
      }

      // if boundary conditions are infinite always set the point
      // if not, check, whether it is inside of domain
      for (unsigned i = 0; i < D; ++i) {
        if (grid.getBoundaryConditions(i) != hrleGrid<D>::INFINITE_BOUNDARY) {
          if (currentIndex[i] > grid.getMaxBounds(i) ||
              currentIndex[i] < grid.getMinBounds(i)) {
            setPoint = false;
          }
        }
      }

      if (setPoint) {
        // Add defined point as it appears in the list
        domain.insertNextDefinedPoint(0, currentIndex, *valueIt);

        // determine signs for next undefined runs
        {
          bool changeSign = false;
          for (int i = D - 1; i >= 0; --i) {
            changeSign = changeSign || (currentIndex[i] > lastIndex[i]);
            if (changeSign) {
              signs[i] = *valueIt < 0;
              lastIndex[i] = currentIndex[i];
            }
          }
        }
      }

      hrleVectorType<hrleIndexType, D> nextIndex;

      ++pointDataIt;
      ++valueIt;

      // choose correct next index
      if (pointDataIt == pointDataEnd) {
        nextIndex = grid.getMaxGridPoint();
        nextIndex[D - 1]++;
      } else {
        for (unsigned i = 0; i < D; ++i) {
          nextIndex[i] = std::round(pointDataIt->at(i) / gridDelta);
        }
      }

      // move current index by one grid spacing and see if the next
      // point has the same index, if not, there must be an undefined
      // run inbetween
      for (int q = 0; q < D; q++) {
        hrleVectorType<hrleIndexType, D> tmp = currentIndex;
        tmp[q]++;
        if (tmp[q] > grid.getMaxGridPoint(q))
          continue;
        for (int r = 0; r < q; ++r)
          tmp[r] = grid.getMinGridPoint(r);

        if (tmp >= nextIndex)
          break;

        domain.insertNextUndefinedPoint(0, tmp,
                                        signs[q] ? lsDomain<T, D>::NEG_VALUE
                                                 : lsDomain<T, D>::POS_VALUE);
      }
    }

    domain.finalize();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsFromMesh)

#endif // LS_FROM_MESH_HPP