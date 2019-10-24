#ifndef LS_MAKE_GEOMETRY_HPP
#define LS_MAKE_GEOMETRY_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleIndexType.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsFromExplicitMesh.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

/// Enumeration for the different types of
/// geometries supported by lsMakeGeometry
enum struct lsMakeGeometryEnum : unsigned {
  SPHERE = 0,
  PLANE = 1,
  BOX = 2
};

/// Create level sets describing basic geometric forms.
template <class T, int D> class lsMakeGeometry {
  typedef typename lsDomain<T, D>::PointValueVectorType pointDataType;

  lsDomain<T, D> *levelSet;
  lsMakeGeometryEnum geometry = lsMakeGeometryEnum::SPHERE;
  const double numericEps = 1e-9;

public:
  lsMakeGeometry(lsDomain<T, D> &passedLevelSet, lsMakeGeometryEnum passedGeometry = lsMakeGeometryEnum::SPHERE) : levelSet(&passedLevelSet), geometry(passedGeometry) {}

  void setLevelSet(lsDomain<T, D> &passedlsDomain) {
    levelSet = &passedlsDomain;
  }

  void setGeometry(lsMakeGeometryEnum passedGeometry) {
    geometry = passedGeometry;
  }

  template <class V> void makeSphere(V origin, T radius, int width = 2) {
    if(levelSet == nullptr) {
      lsMessage::getInstance().addWarning("No level set was passed to lsMakeGeometry.").print();
      return;
    }

    // TODO, this is a stupid algorithm and scales with volume, which is madness
    auto &grid = levelSet->getGrid();
    hrleCoordType gridDelta = grid.getGridDelta();

    hrleVectorType<hrleIndexType, D> index(grid.getMinBounds());
    hrleVectorType<hrleIndexType, D> endIndex(grid.getMaxBounds());

    for (unsigned i = 0; i < D; ++i) {
      if (grid.getBoundaryConditions(i) ==
          hrleGrid<D>::boundaryType::INFINITE_BOUNDARY) {
        index[i] = (origin[i] - radius) / gridDelta - 1;
        endIndex[i] = (origin[i] + radius) / gridDelta + 1;
      }
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

    levelSet->insertPoints(pointData);
    levelSet->getDomain().segment();
    levelSet->finalize(width);
  }

  /// Creates a plane containing the point origin, with
  /// the plane normal given by normal
  template <class V> void makePlane(const V origin, const V passedNormal) {
    if(levelSet == nullptr) {
      lsMessage::getInstance().addWarning("No level set was passed to lsMakeGeometry.").print();
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

    // find minimum and maximum points in infinite direction
    // there are 2*(D-1) points in the corners of the simulation domain
    std::vector<hrleVectorType<T, 3>> cornerPoints;
    cornerPoints.resize(2 * (D - 1));

    // cyclic permutations
    unsigned j = (i + 1) % D;
    unsigned k = (i + 2) % D;

    double minCoord[2];
    double maxCoord[2];
    for (unsigned n = 0; n < D - 1; ++n) {
      minCoord[n] = gridDelta * grid.getMinBounds((i + n + 1) % D);
      maxCoord[n] = gridDelta * grid.getMaxBounds((i + n + 1) % D);
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
      hrleVectorType<unsigned, 2> line(0, 1);
      mesh.insertNextLine(line);
    } else {
      hrleVectorType<unsigned, 3> triangle(0, 1, 2);
      mesh.insertNextTriangle(triangle);
      triangle = hrleVectorType<unsigned, 3>(0, 3, 1);
      mesh.insertNextTriangle(triangle);
    }
    // now convert mesh to levelset
    lsFromExplicitMesh<T, D>(*levelSet, mesh).apply();
  }

  // This function creates a box starting in minCorner spanning to maxCorner
  template <class V> void makeBox(V minCorner, V maxCorner) {
    if(levelSet == nullptr) {
      lsMessage::getInstance().addWarning("No level set was passed to lsMakeGeometry.").print();
      return;
    }
    // draw all triangles for the surface and then import from the mesh
    std::vector<hrleVectorType<T, 3>> corners;
    corners.resize(std::pow(2, D), hrleVectorType<T, 3>(T(0)));

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
      hrleVectorType<unsigned, 2> lines[4];
      lines[0] = hrleVectorType<unsigned, 2>(0, 1);
      lines[1] = hrleVectorType<unsigned, 2>(1, 3);
      lines[2] = hrleVectorType<unsigned, 2>(3, 2);
      lines[3] = hrleVectorType<unsigned, 2>(2, 0);
      for (unsigned i = 0; i < 4; ++i)
        mesh.insertNextLine(lines[i]);
    } else {
      hrleVectorType<unsigned, 3> triangles[12];
      triangles[0] = hrleVectorType<unsigned, 3>(0, 3, 1);
      triangles[1] = hrleVectorType<unsigned, 3>(0, 2, 3);
      triangles[2] = hrleVectorType<unsigned, 3>(0, 1, 5);
      triangles[3] = hrleVectorType<unsigned, 3>(0, 5, 4);
      triangles[4] = hrleVectorType<unsigned, 3>(0, 4, 2);
      triangles[5] = hrleVectorType<unsigned, 3>(4, 6, 2);
      triangles[6] = hrleVectorType<unsigned, 3>(7, 6, 4);
      triangles[7] = hrleVectorType<unsigned, 3>(7, 4, 5);
      triangles[8] = hrleVectorType<unsigned, 3>(7, 2, 6);
      triangles[9] = hrleVectorType<unsigned, 3>(7, 3, 2);
      triangles[10] = hrleVectorType<unsigned, 3>(1, 3, 5);
      triangles[11] = hrleVectorType<unsigned, 3>(3, 7, 5);
      for (unsigned i = 0; i < 12; ++i)
        mesh.insertNextTriangle(triangles[i]);
    }

    // now convert mesh to levelset
    lsFromExplicitMesh<T, D>(*levelSet, mesh).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsMakeGeometry)

#endif // LS_MAKE_GEOMETRY_HPP
