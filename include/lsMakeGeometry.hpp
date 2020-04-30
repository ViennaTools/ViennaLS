#ifndef LS_MAKE_GEOMETRY_HPP
#define LS_MAKE_GEOMETRY_HPP

#include <cassert>

#include <lsPreCompileMacros.hpp>

#include <hrleIndexType.hpp>
#include <hrleVectorType.hpp>

#include <lsConvexHull.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

/// Create level sets describing basic geometric forms.
template <class T, int D> class lsMakeGeometry {
  typedef typename lsDomain<T, D>::PointValueVectorType pointDataType;

  /// Enumeration for the different types of
  /// geometries supported by lsMakeGeometry
  enum struct lsGeometryEnum : unsigned {
    SPHERE = 0,
    PLANE = 1,
    BOX = 2,
    CUSTOM = 3
  };

  lsDomain<T, D> *levelSet;
  lsGeometryEnum geometry = lsGeometryEnum::SPHERE;
  const lsSphere<T, D> *sphere = nullptr;
  const lsPlane<T, D> *plane = nullptr;
  const lsBox<T, D> *box = nullptr;
  lsPointCloud<T, D> *pointCloud = nullptr;
  const double numericEps = 1e-9;
  bool ignoreBoundaryConditions = false;

public:
  lsMakeGeometry() {}

  lsMakeGeometry(lsDomain<T, D> &passedLevelSet) : levelSet(&passedLevelSet) {}

  lsMakeGeometry(lsDomain<T, D> &passedLevelSet,
                 const lsSphere<T, D> &passedSphere)
      : levelSet(&passedLevelSet), sphere(&passedSphere) {
    geometry = lsGeometryEnum::SPHERE;
  }

  lsMakeGeometry(lsDomain<T, D> &passedLevelSet,
                 const lsPlane<T, D> &passedPlane)
      : levelSet(&passedLevelSet), plane(&passedPlane) {
    geometry = lsGeometryEnum::PLANE;
  }

  lsMakeGeometry(lsDomain<T, D> &passedLevelSet, const lsBox<T, D> &passedBox)
      : levelSet(&passedLevelSet), box(&passedBox) {
    geometry = lsGeometryEnum::BOX;
  }

  lsMakeGeometry(lsDomain<T, D> &passedLevelSet,
                 lsPointCloud<T, D> &passedPointCloud)
      : levelSet(&passedLevelSet), pointCloud(&passedPointCloud) {
    geometry = lsGeometryEnum::CUSTOM;
  }

  void setLevelSet(lsDomain<T, D> &passedlsDomain) {
    levelSet = &passedlsDomain;
  }

  /// Set sphere as geometry to be created in the level set.
  void setGeometry(lsSphere<T, D> &passedSphere) {
    sphere = &passedSphere;
    geometry = lsGeometryEnum::SPHERE;
  }

  /// Set a plane to be created in the level set.
  void setGeometry(lsPlane<T, D> &passedPlane) {
    plane = &passedPlane;
    geometry = lsGeometryEnum::PLANE;
  }

  /// Set a box to be created in the level set.
  void setGeometry(lsBox<T, D> &passedBox) {
    box = &passedBox;
    geometry = lsGeometryEnum::BOX;
  }

  /// Set a point cloud, which is used to create
  /// a geometry from its convex hull.
  void setGeometry(lsPointCloud<T, D> &passedPointCloud) {
    pointCloud = &passedPointCloud;
    geometry = lsGeometryEnum::CUSTOM;
  }

  void setIgnoreBoundaryConditions(bool &passedIgnoreBoundaryConditions) {
    ignoreBoundaryConditions = passedIgnoreBoundaryConditions;
  }

  void apply() {
    switch (geometry) {
    case lsGeometryEnum::SPHERE:
      if (sphere == nullptr) {
        lsMessage::getInstance()
            .addWarning("No lsSphere supplied to lsMakeGeometry. Not creating "
                        "geometry.")
            .print();
      }
      makeSphere(sphere->origin, sphere->radius);
      break;
    case lsGeometryEnum::PLANE:
      if (plane == nullptr) {
        lsMessage::getInstance()
            .addWarning(
                "No lsPlane supplied to lsMakeGeometry. Not creating geometry.")
            .print();
      }
      makePlane(plane->origin, plane->normal);
      break;
    case lsGeometryEnum::BOX:
      if (box == nullptr) {
        lsMessage::getInstance()
            .addWarning(
                "No lsBox supplied to lsMakeGeometry. Not creating geometry.")
            .print();
      }
      makeBox(box->minCorner, box->maxCorner);
      break;
    case lsGeometryEnum::CUSTOM:
      if (pointCloud == nullptr) {
        lsMessage::getInstance()
            .addWarning("No lsPointCloud supplied to lsMakeGeometry. Not "
                        "creating geometry.")
            .print();
      }
      makeCustom(pointCloud);
      break;
    default:
      lsMessage::getInstance()
          .addWarning("Invalid geometry type was specified for lsMakeGeometry. "
                      "Not creating geometry.")
          .print();
    }
  }

private:
  void makeSphere(hrleVectorType<T, D> origin, T radius, int width = 2) {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsMakeGeometry.")
          .print();
      return;
    }
    if (sphere == nullptr) {
      lsMessage::getInstance()
          .addWarning("No sphere was passed to lsMakeGeometry.")
          .print();
      return;
    }

    // TODO, this is a stupid algorithm and scales with volume, which is madness
    auto &grid = levelSet->getGrid();
    hrleCoordType gridDelta = grid.getGridDelta();

    // calculate indices from sphere size
    hrleVectorType<hrleIndexType, D> index;
    hrleVectorType<hrleIndexType, D> endIndex;

    for (unsigned i = 0; i < D; ++i) {
      index[i] = (origin[i] - radius) / gridDelta - 1;
      endIndex[i] = (origin[i] + radius) / gridDelta + 1;
    }

    const T valueLimit = width * 0.5 * gridDelta;
    const T radius2 = radius * radius;

    pointDataType pointData;
    const hrleVectorType<hrleIndexType, D> minIndex = index;

    while (index < endIndex) {
      // take shortest manhatten distance to gridline intersection
      T distance = std::numeric_limits<T>::max();
      for (unsigned i = 0; i < D; ++i) {
        T y = (index[(i + 1) % D] * gridDelta) - origin[(i + 1) % D];
        T z = 0;
        if (D == 3)
          z = (index[(i + 2) % D] * gridDelta) - origin[(i + 2) % D];
        T x = radius2 - y * y - z * z;
        if (x < 0.)
          continue;
        T dirRadius =
            std::abs((index[i] * gridDelta) - origin[i]) - std::sqrt(x);
        if (std::abs(dirRadius) < std::abs(distance))
          distance = dirRadius;
      }

      if (std::abs(distance) <= valueLimit + 1e-10) {
        pointData.push_back(std::make_pair(index, distance / gridDelta));
      }
      int dim = 0;
      for (; dim < D - 1; ++dim) {
        if (index[dim] < endIndex[dim])
          break;
        index[dim] = minIndex[dim];
      }
      ++index[dim];
    }

    // Mirror indices correctly into domain, unless boundary conditions
    // are ignored
    if (!ignoreBoundaryConditions) {
      for (unsigned i = 0; i < pointData.size(); ++i) {
        for (unsigned j = 0; j < D; ++j) {
          if (grid.isBoundaryPeriodic(j)) {
            pointData[i].first[j] =
                grid.globalIndex2LocalIndex(j, pointData[i].first[j]);
          }
        }
      }
    }

    levelSet->insertPoints(pointData);
    levelSet->getDomain().segment();
    levelSet->finalize(width);
  }

  /// Creates a plane containing the point origin, with
  /// the plane normal given by normal
  void makePlane(hrleVectorType<T, D> origin,
                 hrleVectorType<T, D> passedNormal) {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsMakeGeometry.")
          .print();
      return;
    }
    if (plane == nullptr) {
      lsMessage::getInstance()
          .addWarning("No plane was passed to lsMakeGeometry.")
          .print();
      return;
    }

    auto &grid = levelSet->getGrid();
    hrleCoordType gridDelta = grid.getGridDelta();

    // normalise passedNormal
    double modulus = 0.;
    hrleVectorType<double, D> normal(passedNormal);
    for (unsigned i = 0; i < D; ++i) {
      modulus += normal[i] * normal[i];
    }
    modulus = std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      normal[i] /= modulus;
    }

    // check that boundary conditions are correct
    unsigned i = 0;
    bool infDimSet = false;
    for (unsigned n = 0; n < D; ++n) {
      if (grid.getBoundaryConditions(n) ==
          hrleGrid<D>::boundaryType::INFINITE_BOUNDARY) {
        if (!infDimSet) {
          i = n;
          infDimSet = true;
        } else {
          lsMessage::getInstance().addError(
              "Planes can only be created with one Infinite Boundary "
              "Condition. More than one found!");
        }
      }
    }
    if (!infDimSet) {
      lsMessage::getInstance().addError("Planes require exactly one Infinite "
                                        "Boundary Condition. None found!");
    }

    if (passedNormal[i] == 0.) {
      lsMessage::getInstance().addError(
          "lsMakeGeometry: Plane cannot be parallel to Infinite Boundary "
          "direction!");
    }

    // find minimum and maximum points in infinite direction
    // there are 2*(D-1) points in the corners of the simulation domain
    std::vector<std::array<double, 3>> cornerPoints;
    cornerPoints.resize(2 * (D - 1));

    // cyclic permutations
    unsigned j = (i + 1) % D;
    unsigned k = (i + 2) % D;

    double minCoord[2];
    double maxCoord[2];
    // Find grid boundaries, there used to be a +-1 for the coords.
    // If an error pops up here, probably has to do with that.
    // But if +-1 is added here, the boundaries are exceeded and
    // the correct boundary conditions will add stray points for
    // tilted planes in lsFromSurfaceMesh later on.
    for (unsigned n = 0; n < D - 1; ++n) {
      minCoord[n] = gridDelta * (grid.getMinIndex((i + n + 1) % D));
      maxCoord[n] = gridDelta * (grid.getMaxIndex((i + n + 1) % D));
    }

    // set corner points
    cornerPoints[0][j] = minCoord[0];
    cornerPoints[1][j] = maxCoord[0];

    if (D == 3) {
      cornerPoints[0][k] = minCoord[1];
      cornerPoints[1][k] = maxCoord[1];

      cornerPoints[2][j] = minCoord[0];
      cornerPoints[2][k] = maxCoord[1];
      cornerPoints[3][j] = maxCoord[0];
      cornerPoints[3][k] = minCoord[1];
    }

    // now find i coordinate of points
    lsMesh mesh;

    for (unsigned n = 0; n < cornerPoints.size(); ++n) {
      double numerator = (cornerPoints[n][j] - origin[j]) * normal[j];
      if (D == 3)
        numerator += (cornerPoints[n][k] - origin[k]) * normal[k];
      else
        cornerPoints[n][2] = 0.;
      cornerPoints[n][i] = origin[i] - numerator / normal[i];
      mesh.insertNextNode(cornerPoints[n]);
    }

    if (D == 2) {
      std::array<unsigned, 2> line = {0, 1};
      if (normal[i] < 0.)
        std::swap(line[0], line[1]);
      mesh.insertNextLine(line);
    } else {
      std::array<unsigned, 3> triangle = {0, 1, 2};
      if (normal[i] < 0.)
        std::swap(triangle[0], triangle[1]);
      mesh.insertNextTriangle(triangle);
      triangle = {0, 3, 1};
      if (normal[i] < 0.)
        std::swap(triangle[0], triangle[1]);
      mesh.insertNextTriangle(triangle);
    }
    // now convert mesh to levelset
    lsFromSurfaceMesh<T, D>(*levelSet, mesh).apply();
  }

  // This function creates a box starting in minCorner spanning to maxCorner
  void makeBox(hrleVectorType<T, D> minCorner, hrleVectorType<T, D> maxCorner) {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsMakeGeometry.")
          .print();
      return;
    }
    if (box == nullptr) {
      lsMessage::getInstance()
          .addWarning("No box was passed to lsMakeGeometry.")
          .print();
      return;
    }

    // draw all triangles for the surface and then import from the mesh
    std::vector<std::array<double, 3>> corners;
    corners.resize(std::pow(2, D), {0, 0, 0});

    // first corner is the minCorner
    for (unsigned i = 0; i < D; ++i)
      corners[0][i] = minCorner[i];

    // last corner is maxCorner
    for (unsigned i = 0; i < D; ++i)
      corners.back()[i] = maxCorner[i];

    // calculate all missing corners
    corners[1] = corners[0];
    corners[1][0] = corners.back()[0];

    corners[2] = corners[0];
    corners[2][1] = corners.back()[1];

    if (D == 3) {
      corners[3] = corners.back();
      corners[3][2] = corners[0][2];

      corners[4] = corners[0];
      corners[4][2] = corners.back()[2];

      corners[5] = corners.back();
      corners[5][1] = corners[0][1];

      corners[6] = corners.back();
      corners[6][0] = corners[0][0];
    }

    // now add all corners to mesh
    lsMesh mesh;
    for (unsigned i = 0; i < corners.size(); ++i) {
      mesh.insertNextNode(corners[i]);
    }

    if (D == 2) {
      std::array<unsigned, 2> lines[4] = {{0, 2}, {2, 3}, {3, 1}, {1, 0}};
      for (unsigned i = 0; i < 4; ++i)
        mesh.insertNextLine(lines[i]);
    } else {
      std::array<unsigned, 3> triangles[12] = {
          {0, 3, 1}, {0, 2, 3}, {0, 1, 5}, {0, 5, 4}, {0, 4, 2}, {4, 6, 2},
          {7, 6, 4}, {7, 4, 5}, {7, 2, 6}, {7, 3, 2}, {1, 3, 5}, {3, 7, 5}};
      for (unsigned i = 0; i < 12; ++i)
        mesh.insertNextTriangle(triangles[i]);
    }

    // now convert mesh to levelset
    lsFromSurfaceMesh<T, D>(*levelSet, mesh, ignoreBoundaryConditions).apply();
  }

  void makeCustom(lsPointCloud<T, D> *pointCloud) {
    // create mesh from point cloud
    lsMesh mesh;
    lsConvexHull<T, D>(mesh, *pointCloud).apply();

    // read mesh from surface
    lsFromSurfaceMesh<T, D>(*levelSet, mesh, ignoreBoundaryConditions).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsMakeGeometry)

#endif // LS_MAKE_GEOMETRY_HPP
