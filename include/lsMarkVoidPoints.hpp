#ifndef LS_MARK_VOID_POINTS_HPP
#define LS_MARK_VOID_POINTS_HPP

#include <hrleSparseStarIterator.hpp>

#include <lsPreCompileMacros.hpp>

#include <lsDomain.hpp>
#include <lsGraph.hpp>

/// Enumeration describing which connected component to use
/// as top surface during void point detection.
/// All others points will be set as void poitns.
enum struct lsVoidTopSurfaceEnum : unsigned {
  LEX_LOWEST = 0,
  LEX_HIGHEST = 1,
  LARGEST = 2,
  SMALLEST = 3
};

/// This class is used to mark points of the level set
/// which are enclosed in a void.
template <class T, int D> class lsMarkVoidPoints {
  using IndexType = std::size_t;

  lsSmartPointer<lsDomain<T, D>> domain = nullptr;
  bool reverseVoidDetection = false;
  bool saveComponents = false;
  bool detectLargestSurface = false;

  // two points are connected if they have the same sign
  bool areConnected(const T &value1, const T &value2) {
    return (value1 >= 0) == (value2 >= 0);
  }

  std::vector<IndexType> mergeComponentCounts(const std::vector<IndexType>& components, const std::vector<IndexType>& pointsPerComponent) {
    // find number of connected components after merge
    // TODO: the last element in the components vector is very likely
    // the highest index, so could just set it instead of checking all values
    unsigned highestIndex = 0;
    for(auto it : components) {
      if(highestIndex < it) {
        highestIndex = it;
      }
    }

    // merge pointsPerComponent acording to the connected components
    std::vector<IndexType> pointsPerConnected;
    pointsPerConnected.resize(highestIndex + 1, 0);
    for(unsigned i = 0; i < components.size(); ++i) {
      pointsPerConnected[components[i]] += pointsPerComponent[i];
    }

    return pointsPerConnected;
  }

  IndexType calculateTopID(const std::vector<IndexType>& components, const std::vector<IndexType>& pointsPerConnected) {
    // check which component has the most points
    IndexType topId = 0;
    // use first component, which contains more than 0 points
    while(topId < pointsPerConnected.size() && pointsPerConnected[topId] == 0) {
      ++topId;
    }
    for(unsigned i = topId + 1; i < pointsPerConnected.size(); ++i) {
      if((pointsPerConnected[topId] < pointsPerConnected[i]) != reverseVoidDetection) {
        topId = i;
      }
    }

    return topId;
  }

public:
  lsMarkVoidPoints() {}

  lsMarkVoidPoints(lsSmartPointer<lsDomain<T, D>> passedlsDomain,
                   bool passedReverseVoidDetection = false)
      : domain(passedlsDomain),
        reverseVoidDetection(passedReverseVoidDetection) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    domain = passedlsDomain;
  }

  /// Set whether the "top" level set should be the most positive(default)
  /// connected chain of level set values, or the most negative.
  /// Most positive/negative refers to the lexicographical ordering
  /// of the coordinate of the point.
  void setReverseVoidDetection(bool passedReverseVoidDetection) {
    reverseVoidDetection = passedReverseVoidDetection;
  }

  /// Set whether the number of points of one connected surface
  /// should be used to detect void points. Defaults to false.
  /// If this is set to true, the largest connected surface will be 
  /// kept and all other grid points marked as void points.
  /// By setting reverseVoidDetection to true, the smallest
  /// surface will be used instead.
  void setDetectLargestSurface(bool passedDetect) {
    detectLargestSurface = passedDetect;
  }

  /// Set which connected component to use as the top surface
  /// and mark all other components as void points.
  void setVoidTopSurface(lsVoidTopSurfaceEnum topSurface) {

    switch (topSurface)
    {
      case lsVoidTopSurfaceEnum::LEX_LOWEST:
        reverseVoidDetection = true;
        detectLargestSurface = false;
        break;
      case lsVoidTopSurfaceEnum::LEX_HIGHEST:
        reverseVoidDetection = false;
        detectLargestSurface = false;
        break;
      case lsVoidTopSurfaceEnum::LARGEST:
        reverseVoidDetection = false;
        detectLargestSurface = true;
        break;
      case lsVoidTopSurfaceEnum::SMALLEST:
        reverseVoidDetection = true;
        detectLargestSurface = true;
        break;
    
      default:
        lsMessage::getInstance().addWarning("lsMarkVoidPoints: Invalid lsVoidTopSurfaceEnum set. Using default values.").print();
        reverseVoidDetection = false;
        detectLargestSurface = false;
        break;
    }
  }

  /// Set whether the connected component IDs used to generate the void
  /// points should be saved. Ech point is assigned a component ID
  /// denoting which other points it is connected to.
  void setSaveComponentIds(bool scid) { saveComponents = scid; }

  void apply() {
    lsInternal::lsGraph graph;

    std::vector<std::vector<std::vector<IndexType>>> componentList(
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
    std::vector<IndexType> pointsPerComponent;

    // cycle through and set up the graph to get connectivity information
    for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
             neighborIt(domain->getDomain());
         !neighborIt.isFinished(); neighborIt.next()) {
      auto &center = neighborIt.getCenter();

      // component id of the current run
      auto &currentComponentId =
          componentList[center.getSegmentId()][center.getLevel()]
                       [center.getRunTypePosition()];

      // -1 means it is not set yet
      if (currentComponentId == -1) {
        for (int k = 0; k < 2 * D; ++k) {
          auto &neighbor = neighborIt.getNeighbor(k);
          const auto neighborComponentId =
              componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                           [neighbor.getRunTypePosition()];
          // if neighbor is already defined
          // set current component to neighbor component
          // if they are connected
          if (neighborComponentId != -1) {
            if (areConnected(center.getValue(), neighbor.getValue())) {
              currentComponentId = neighborComponentId;
              if(center.getValue() >= 0.)
                ++pointsPerComponent[currentComponentId];
              break;
            }
          }
        }
      }

      // it is still not set, so add new vertex
      if (currentComponentId == -1) {
        currentComponentId = numberOfComponents;
        pointsPerComponent.push_back((center.getValue() > 0.)?1:0);
        graph.insertNextVertex();
        ++numberOfComponents;
      }

      // check if edge can be set
      for (int k = 0; k < 2 * D; ++k) {
        auto &neighbor = neighborIt.getNeighbor(k);
        auto &neighborComponentId =
            componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                         [neighbor.getRunTypePosition()];
        if (areConnected(center.getValue(), neighbor.getValue())) {
          // if neighbor is already set
          if (neighborComponentId != -1) {
            // if neighbor is part of different component
            if (currentComponentId != neighborComponentId) {
              graph.insertNextEdge(currentComponentId, neighborComponentId);
            }
          } else {
            neighborComponentId = currentComponentId;
            if(neighbor.getValue() >= 0.)
              ++pointsPerComponent[neighborComponentId];
          }
        }
      }
    }

    auto components = graph.getConnectedComponents();

    int topComponent = (reverseVoidDetection) ? components[0] : components.back();

    // identify which layer to keep
    {
      auto pointsPerConnected = mergeComponentCounts(components, pointsPerComponent);

      // find largest connected surface
      if(detectLargestSurface) {
        topComponent = calculateTopID(components, pointsPerConnected);
      } else { // if component does not contain points, take the next one
        while(pointsPerConnected[topComponent] == 0) {
          if(reverseVoidDetection)
            ++topComponent;
          else
            --topComponent;
        }
      }
    }

    std::vector<T> voidPointMarkers;
    voidPointMarkers.resize(domain->getNumberOfPoints());

    std::vector<T> componentMarkers;
    if (saveComponents)
      componentMarkers.resize(domain->getNumberOfPoints());

    // cycle through again to set correct voidPointMarkers
    for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
             neighborIt(domain->getDomain());
         !neighborIt.isFinished(); neighborIt.next()) {
      auto center = neighborIt.getCenter();

      if (!center.isDefined())
        continue;

      // if it is positive, just check if it is part of the top component
      if (center.getValue() >= 0) {
        const int &oldComponentId = componentList[center.getSegmentId()][0]
                                                 [center.getRunTypePosition()];
        voidPointMarkers[center.getPointId()] =
            (components[oldComponentId] != topComponent);
      } else {
        // if it is negative, check all neighbours (with different sign),
        // because they are part of of a positive component, which might
        // be the top
        unsigned k;
        for (k = 0; k < 2 * D; ++k) {
          auto &neighbor = neighborIt.getNeighbor(k);
          if(std::signbit(neighbor.getValue()) == std::signbit(center.getValue()))
            continue;
          const int &oldneighborComponentId =
              componentList[neighbor.getSegmentId()][neighbor.getLevel()]
                           [neighbor.getRunTypePosition()];
          if (components[oldneighborComponentId] == topComponent) {
            break;
          }
        }
        voidPointMarkers[center.getPointId()] = (k == 2 * D);
      }

      if (saveComponents) {
        componentMarkers[center.getPointId()] =
            componentList[center.getSegmentId()][0]
                         [center.getRunTypePosition()];
      }
    }

    auto &pointData = domain->getPointData();
    auto voidMarkersPointer = pointData.getScalarData("VoidPointMarkers");
    // if vector data does not exist
    if (voidMarkersPointer == nullptr) {
      pointData.insertNextScalarData(voidPointMarkers, "VoidPointMarkers");
    } else {
      *voidMarkersPointer = std::move(voidPointMarkers);
    }

    if (saveComponents) {
      auto componentMarkersPointer =
          pointData.getScalarData("ConnectedComponentId");
      // if vector data does not exist
      if (componentMarkersPointer == nullptr) {
        pointData.insertNextScalarData(componentMarkers,
                                       "ConnectedComponentId");
      } else {
        *componentMarkersPointer = std::move(componentMarkers);
      }
    }
  }
};

#endif // LS_MARK_VOID_POINTS_HPP
