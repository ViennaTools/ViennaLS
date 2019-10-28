#ifndef LS_MARK_VOID_POINTS_HPP
#define LS_MARK_VOID_POINTS_HPP

#include <hrleSparseStarIterator.hpp>

#include <lsPreCompileMacros.hpp>

#include <lsDomain.hpp>
#include <lsGraph.hpp>

/// This class is used to mark points of the level set
/// which are enclosed in a void.
template <class T, int D> class lsMarkVoidPoints {
  lsDomain<T, D> *domain = nullptr;
  bool reverseVoidDetection = false;

  // two points are connected if they have the same sign
  bool areConnected(const T &value1, const T &value2) {
    return (value1 >= 0) == (value2 >= 0);
  }

public:
  lsMarkVoidPoints(lsDomain<T, D> &passedlsDomain,
                   bool passedReverseVoidDetection = false)
      : domain(&passedlsDomain),
        reverseVoidDetection(passedReverseVoidDetection) {}

  void setLevelSet(lsDomain<T, D> &passedlsDomain) { domain = &passedlsDomain; }

  void setReverseVoidDetection(bool passedReverseVoidDetection) {
    reverseVoidDetection = passedReverseVoidDetection;
  }

  void apply() {
    lsInternal::lsGraph graph;

    // vector<int> with [segmentId][dimension+1][runtype] to save connectiviies
    std::vector<std::vector<std::vector<int>>> componentList(
        domain->getNumberOfSegments());
    for (unsigned segmentId = 0; segmentId < domain->getNumberOfSegments();
         ++segmentId) {
      componentList[segmentId].resize(D + 1);
      for (int dim = -1; dim < D; ++dim) {
        componentList[segmentId][dim + 1].resize(
            domain->getDomain().getNumberOfRuns(segmentId, dim), -1);
      }
    }

    std::size_t numberOfComponents = 0;

    // cycle through and set up the graph to get connectivity information
    for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
             neighborIt(domain->getDomain());
         !neighborIt.isFinished(); neighborIt.next()) {
      auto &center = neighborIt.getCenter();

      // component id of the current run
      // std::cout << center.getSegmentId() << ", " << center.getLevel() << ", "
      // << center.getRunTypePosition() << std::endl;
      int &currentComponentId =
          componentList[center.getSegmentId()][center.getLevel()]
                       [center.getRunTypePosition()];

      // -1 means it is not set yet
      if (currentComponentId == -1) {
        for (int k = 0; k < 2 * D; ++k) {
          auto &neighbor = neighborIt.getNeighbor(k);
          const int &neighborComponentId =
              componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                           [neighbor.getRunTypePosition()];
          // if neighbor is already defined
          // set current component to neighbor component
          // if they are connected
          if (neighborComponentId != -1) {
            if (areConnected(center.getValue(), neighbor.getValue())) {
              currentComponentId = neighborComponentId;
              break;
            }
          }
        }
      }

      // it is still not set, so add new vertex
      if (currentComponentId == -1) {
        currentComponentId = numberOfComponents;
        // std::cout << "Inserting vertex: " << currentComponentId << std::endl;
        graph.insertNextVertex();
        ++numberOfComponents;
      }

      // check if edge can be set
      for (int k = 0; k < 2 * D; ++k) {
        auto &neighbor = neighborIt.getNeighbor(k);
        int &neighborComponentId =
            componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                         [neighbor.getRunTypePosition()];
        if (areConnected(center.getValue(), neighbor.getValue())) {
          // if neighbor is already set
          if (neighborComponentId != -1) {
            // if neighbor is part of different component
            if (currentComponentId != neighborComponentId) {
              // std::cout << "Inserting: " << currentComponentId << ", " <<
              // neighborComponentId << std::endl; std::cout << k << ", " <<
              // neighbor.getSegmentId() << ", " << neighbor.getLevel() << ", "
              // << neighbor.getRunTypePosition() << std::endl;
              graph.insertNextEdge(currentComponentId, neighborComponentId);
            }
          } else {
            // std::cout << "Setting neighbor to: " << currentComponentId <<
            // std::endl;
            neighborComponentId = currentComponentId;
          }
        }
      }
    }

    auto components = graph.getConnectedComponents();

    // now need to decide which component to mark as valid
    // always take the component in the most positive direction
    // unless reverse is specified, then take the most negative
    int topComponent =
        (reverseVoidDetection) ? components[0] : components.back();

    std::vector<bool> &voidPointMarkers = domain->getVoidPointMarkers();
    voidPointMarkers.reserve(domain->getNumberOfPoints());

    // cycle through again to set correct voidPointMarkers
    for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
             neighborIt(domain->getDomain());
         !neighborIt.isFinished(); neighborIt.next()) {
      auto center = neighborIt.getCenter();

      if (!center.isDefined())
        continue;

      // TODO: currently this will break if the iterator hits a negative
      // run first
      if (center.getValue() >= 0) {
        const int &oldComponentId = componentList[center.getSegmentId()][0]
                                                 [center.getRunTypePosition()];
        voidPointMarkers.push_back(components[oldComponentId] == topComponent);
      } else {
        unsigned k;
        for (k = 0; k < 2 * D; ++k) {
          auto &neighbor = neighborIt.getNeighbor(k);
          const int &oldneighborComponentId =
              componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                           [neighbor.getRunTypePosition()];
          if (components[oldneighborComponentId] == topComponent) {
            break;
          }
        }
        voidPointMarkers.push_back(k != 2 * D);
      }
    }
  }
};

#endif // LS_MARK_VOID_POINTS_HPP
