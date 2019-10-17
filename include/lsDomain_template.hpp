#ifndef LS_DOMAIN_TEMPLATE_HPP
#define LS_DOMAIN_TEMPLATE_HPP

#include <limits>

#include <hrleDomain.hpp>
#include <hrleFillDomainWithSignedDistance.hpp>
#include <hrleVectorType.hpp>

template <class T, int D> class lsDomain {
public:
  // TYPEDEFS
  typedef T ValueType;
  typedef hrleGrid<D> GridType;
  typedef hrleDomain<T, D> DomainType;
  typedef typename GridType::boundaryType BoundaryType;
  typedef typename std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>>
      PointValueVectorType;

private:
  // PRIVATE MEMBER VARIABLES
  GridType grid;
  DomainType domain;
  int levelSetWidth = 1;
  unsigned numberOfActivePoints;
  std::vector<hrleSizeType> activePointIds;

public:
  // STATIC CONSTANTS
  static constexpr int dimensions = D;

  // PUBLIC MEMBER VARIABLES
  static constexpr T POS_VALUE = std::numeric_limits<T>::max();
  static constexpr T NEG_VALUE = std::numeric_limits<T>::lowest();

  /// initalise an empty infinite lsDomain
  lsDomain(hrleCoordType gridDelta = 1.0) {
    hrleIndexType gridMin[D], gridMax[D];
    BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D; ++i) {
      gridMin[i] = 0;
      gridMax[i] = 0;
      boundaryCons[i] = BoundaryType::INFINITE_BOUNDARY;
    }

    grid = GridType(gridMin, gridMax, gridDelta, boundaryCons);
    domain.deepCopy(grid, DomainType(grid, T(NEG_VALUE)));
  }

  lsDomain(hrleCoordType *bounds, BoundaryType *boundaryConditions,
           hrleCoordType gridDelta = 1.0) {
    hrleIndexType gridMin[D], gridMax[D];
    for (unsigned i = 0; i < D; ++i) {
      gridMin[i] = std::floor(bounds[2 * i] / gridDelta);
      gridMax[i] = std::ceil(bounds[2 * i + 1] / gridDelta);
    }

    grid = GridType(gridMin, gridMax, gridDelta, boundaryConditions);
    domain.deepCopy(grid, DomainType(grid, T(NEG_VALUE)));
  }

  /// initialise lsDomain with domain size "bounds", filled with point/value
  /// pairs in pointData
  lsDomain(PointValueVectorType pointData, hrleCoordType *bounds,
           BoundaryType *boundaryConditions, hrleCoordType gridDelta = 1.0) {
    this->deepCopy(lsDomain(bounds, boundaryConditions, gridDelta));
    hrleFillDomainWithSignedDistance(domain, pointData, T(NEG_VALUE),
                                     T(POS_VALUE));
  }

  lsDomain(GridType passedGrid) : grid(passedGrid) {
    domain.deepCopy(grid, DomainType(grid, T(NEG_VALUE)));
  }

  lsDomain(lsDomain &passedlsDomain) { deepCopy(passedlsDomain); }

  /// this function sets a new levelset width and finalizes the levelset, so it
  /// is ready for use by other algorithms
  void finalize(int newWidth) { levelSetWidth = newWidth; }

  /// this function finalizes the levelset, so it is ready for use by other
  /// algorithms
  void finalize() {}

  /// copy all values of "passedlsDomain" to this lsDomain
  void deepCopy(const lsDomain<T, D> &passedlsDomain) {
    grid = passedlsDomain.grid;
    domain.deepCopy(grid, passedlsDomain.domain);
    levelSetWidth = passedlsDomain.levelSetWidth;
  }

  /// re-initalise lsDomain with the point/value pairs in pointData
  void insertPoints(PointValueVectorType pointData) {
    hrleFillDomainWithSignedDistance(domain, pointData, T(NEG_VALUE),
                                     T(POS_VALUE));
  }

  /// constructs the activePointIds vector which holds the active point ID
  /// of a point at the index PointId
  void calculateActivePointIds() {
    activePointIds.clear();
    hrleSparseIterator<hrleDomain<T, D>> it(domain);
    hrleSizeType currentActiveID = 0;
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      activePointIds.push_back((std::abs(it.getValue()) > 0.5)
                                   ? std::numeric_limits<hrleSizeType>::max()
                                   : currentActiveID++);
    }
    numberOfActivePoints = currentActiveID - 1;
  }

  /// get the activePointId of the point with pointId
  hrleSizeType getActivePointId(hrleSizeType pointId) const {
    return activePointIds[pointId];
  }

  /// get the pointId from the active point with activePointId
  hrleSizeType getPointId(hrleSizeType activePointId) const {
    return std::distance(
        activePointIds.begin(),
        std::find(activePointIds.begin(), activePointIds.end(), activePointId));
  }

  /// get reference to the grid on which the levelset is defined
  const GridType &getGrid() const { return grid; }

  /// get const reference to the underlying hrleDomain data structure
  DomainType &getDomain() { return domain; }

  const DomainType &getDomain() const { return domain; }

  /// returns the number of segments, the levelset is split into.
  /// This is useful for algorithm parallelisation
  unsigned getNumberOfSegments() const { return domain.getNumberOfSegments(); }

  /// returns the number of defined points
  unsigned getNumberOfPoints() const { return domain.getNumberOfPoints(); }

  /// returns the number of defined point with value <= 0.5
  unsigned getNumberOfActivePoints() const { return numberOfActivePoints; }

  int getLevelSetWidth() const { return levelSetWidth; }

  void setLevelSetWidth(int width) { levelSetWidth = width; }

  /// prints basic information and all memebers of the levelset structure
  void print() {
    std::cout << "Grid pointer: " << &grid << std::endl;
    std::cout << "Domain: " << &domain << std::endl;
    std::cout << "DomainSegments: " << std::endl;
    for (unsigned i = 0; i < getNumberOfSegments(); ++i) {
      std::cout << &(domain.getDomainSegment(i)) << std::endl;
    }
    domain.print();
  }
};

#endif // LS_DOMAIN_TEMPLATE_HPP
