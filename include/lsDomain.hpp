#ifndef LS_DOMAIN_HPP
#define LS_DOMAIN_HPP

#include <lsPreCompileMacros.hpp>

#include <limits>

#include <hrleDomain.hpp>
#include <hrleFillDomainWithSignedDistance.hpp>
#include <hrleVectorType.hpp>

#include <lsPointData.hpp>
#include <lsSmartPointer.hpp>

#define LS_DOMAIN_SERIALIZATION_VERSION 0

///  Class containing all information about the level set, including
///  the dimensions of the domain, boundary conditions and all data.
template <class T, int D> class lsDomain {
public:
  // TYPEDEFS
  typedef T ValueType;
  typedef hrleGrid<D> GridType;
  typedef hrleDomain<T, D> DomainType;
  typedef typename GridType::boundaryType BoundaryType;
  typedef typename std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>>
      PointValueVectorType;
  typedef typename std::vector<std::array<T, D>> NormalVectorType;
  typedef lsPointData PointDataType;
  typedef typename std::vector<bool> VoidPointMarkersType;

private:
  // PRIVATE MEMBER VARIABLES
  GridType grid;
  DomainType domain;
  int levelSetWidth = 1;
  PointDataType pointData;
  VoidPointMarkersType voidPointMarkers;

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
    domain.deepCopy(grid, DomainType(grid, T(POS_VALUE)));
  }

  lsDomain(hrleCoordType *bounds, BoundaryType *boundaryConditions,
           hrleCoordType gridDelta = 1.0) {
    hrleIndexType gridMin[D], gridMax[D];
    for (unsigned i = 0; i < D; ++i) {
      gridMin[i] = std::floor(bounds[2 * i] / gridDelta);
      gridMax[i] = std::ceil(bounds[2 * i + 1] / gridDelta);
    }

    grid = GridType(gridMin, gridMax, gridDelta, boundaryConditions);
    domain.deepCopy(grid, DomainType(grid, T(POS_VALUE)));
  }

  lsDomain(std::vector<hrleCoordType> bounds,
           std::vector<unsigned> boundaryConditions,
           hrleCoordType gridDelta = 1.0) {
    BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D; ++i) {
      boundaryCons[i] = static_cast<BoundaryType>(boundaryConditions[i]);
    }
    auto newDomain = lsSmartPointer<lsDomain<T, D>>::New(
        bounds.data(), boundaryCons, gridDelta);
    this->deepCopy(newDomain);
  }

  /// initialise lsDomain with domain size "bounds", filled with point/value
  /// pairs in pointData
  lsDomain(PointValueVectorType pointData, hrleCoordType *bounds,
           BoundaryType *boundaryConditions, hrleCoordType gridDelta = 1.0) {
    auto newDomain = lsSmartPointer<lsDomain<T, D>>::New(
        bounds, boundaryConditions, gridDelta);
    this->deepCopy(newDomain);
    hrleFillDomainWithSignedDistance(domain, pointData, T(NEG_VALUE),
                                     T(POS_VALUE));
  }

  lsDomain(GridType passedGrid) : grid(passedGrid) {
    domain.deepCopy(grid, DomainType(grid, T(POS_VALUE)));
  }

  lsDomain(lsSmartPointer<lsDomain> passedDomain) { deepCopy(passedDomain); }

  /// this function sets a new levelset width and finalizes the levelset, so it
  /// is ready for use by other algorithms
  void finalize(int newWidth) { levelSetWidth = newWidth; }

  /// this function finalizes the levelset, so it is ready for use by other
  /// algorithms
  void finalize() {}

  /// copy all values of "passedlsDomain" to this lsDomain
  void deepCopy(const lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    grid = passedlsDomain->grid;
    domain.deepCopy(grid, passedlsDomain->domain);
    levelSetWidth = passedlsDomain->levelSetWidth;
    pointData = passedlsDomain->pointData;
  }

  /// re-initalise lsDomain with the point/value pairs in pointData
  /// This is similar to lsFromMesh with the difference that pointData
  /// contains (INDEX, Value) pairs, while lsFromMesh expects coordinates
  /// rather than indices
  void insertPoints(PointValueVectorType pointData, bool sort = true) {
    hrleFillDomainWithSignedDistance(domain, pointData, T(NEG_VALUE),
                                     T(POS_VALUE), sort);
  }

  /// get reference to the grid on which the levelset is defined
  const GridType &getGrid() const { return grid; }

  /// get mutable reference to the grid on which the level set is defined
  GridType &getGrid() { return grid; }

  /// get const reference to the underlying hrleDomain data structure
  DomainType &getDomain() { return domain; }

  const DomainType &getDomain() const { return domain; }

  /// returns the number of segments, the levelset is split into.
  /// This is useful for algorithm parallelisation
  unsigned getNumberOfSegments() const { return domain.getNumberOfSegments(); }

  /// returns the number of defined points
  unsigned getNumberOfPoints() const { return domain.getNumberOfPoints(); }

  int getLevelSetWidth() const { return levelSetWidth; }

  void setLevelSetWidth(int width) { levelSetWidth = width; }

  // clear all additional data
  void clearMetaData() { pointData.clear(); }

  /// get reference to point data saved in the level set
  PointDataType &getPointData() { return pointData; }

  const PointDataType &getPointData() const { return pointData; }

  /// get reference to the voidPoints markers for all points
  VoidPointMarkersType &getVoidPointMarkers() { return voidPointMarkers; }

  const VoidPointMarkersType &getVoidPointMarkers() const {
    return voidPointMarkers;
  }

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

  /// Serializes the lsDomain into a binary stream
  std::ostream &serialize(std::ostream &stream) {
    // Save header to identify lsDomain
    stream << "lsDomain";

    // now write format version number
    char formatVersion = LS_DOMAIN_SERIALIZATION_VERSION;
    stream.write(&formatVersion, 1);

    // serialize grid
    grid.serialize(stream);

    // serialize hrleDomain which saves LS values
    domain.serialize(stream);

    // serialize lsDomain members
    // level set width as 32bit uint
    const uint32_t width = levelSetWidth;
    stream.write(reinterpret_cast<const char *>(&width), sizeof(uint32_t));

    // serialize pointData if there is any point data associated with this
    // lsDomain mark whether there is point data or not (1/0)
    char hasPointData = (pointData.empty()) ? 0 : 1;
    stream.write(&hasPointData, 1);
    if (hasPointData == 1) {
      pointData.serialize(stream);
    }

    return stream;
  }

  /// Deserialize lsDomain from binary stream
  std::istream &deserialize(std::istream &stream) {
    // Check identifier
    char identifier[8];
    stream.read(identifier, 8);
    if (std::string(identifier).compare(0, 8, "lsDomain")) {
      lsMessage::getInstance()
          .addWarning(
              "Reading lsDomain from stream failed. Header could not be found.")
          .print();
      return stream;
    }

    // check format version for compatibility
    char formatVersion;
    stream.read(&formatVersion, 1);
    if (formatVersion > LS_DOMAIN_SERIALIZATION_VERSION) {
      lsMessage::getInstance()
          .addWarning(
              "Reading lsDomain of version " + std::to_string(formatVersion) +
              " with reader of version " +
              std::to_string(LS_DOMAIN_SERIALIZATION_VERSION) + " failed.")
          .print();
      return stream;
    }

    // read in the grid
    grid.deserialize(stream);

    // read in the hrleDoamin
    // grid pointer in hrleDomain is implicitly set
    // to the correct grid already since they were
    // initialized together
    domain.deserialize(stream);

    // read in the level set width
    uint32_t width;
    stream.read(reinterpret_cast<char *>(&width), sizeof(uint32_t));
    levelSetWidth = width;

    // check wether there is point data to read
    char hasPointData;
    stream.read(&hasPointData, 1);
    if (hasPointData == 1) {
      pointData.clear();
      pointData.deserialize(stream);
    }

    return stream;
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsDomain)

#endif // LS_DOMAIN_HPP
