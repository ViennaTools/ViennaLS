#pragma once

#include <cassert>

#include <lsPreCompileMacros.hpp>

#include <hrleTypes.hpp>

#include <lsConvexHull.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMesh.hpp>
#include <lsTransformMesh.hpp>

#include <vcVectorType.hpp>

#ifndef NDEBUG
#include <lsVTKWriter.hpp>
#endif

namespace viennals {

using namespace viennacore;

/// Create level sets describing basic geometric forms.
template <class T, int D> class MakeGeometry {
  typedef typename Domain<T, D>::PointValueVectorType pointDataType;

  /// Enumeration for the different types of
  /// geometries supported by MakeGeometry
  enum struct GeometryEnum : unsigned {
    SPHERE = 0,
    PLANE = 1,
    BOX = 2,
    CUSTOM = 3,
    CYLINDER = 4
  };

  SmartPointer<Domain<T, D>> levelSet;
  GeometryEnum geometry = GeometryEnum::SPHERE;
  SmartPointer<Sphere<T, D>> sphere;
  SmartPointer<Plane<T, D>> plane;
  SmartPointer<Box<T, D>> box;
  SmartPointer<Cylinder<T, D>> cylinder;
  SmartPointer<PointCloud<T, D>> pointCloud;
  const double numericEps = 1e-9;
  // bool ignoreBoundaryConditions = false;
  std::array<bool, 3> ignoreBoundaryConditions{false, false, false};

public:
  MakeGeometry() = default;

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet,
               SmartPointer<Sphere<T, D>> passedSphere)
      : levelSet(passedLevelSet), sphere(passedSphere) {
    geometry = GeometryEnum::SPHERE;
  }

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet,
               SmartPointer<Plane<T, D>> passedPlane)
      : levelSet(passedLevelSet), plane(passedPlane) {
    geometry = GeometryEnum::PLANE;
  }

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet,
               SmartPointer<Box<T, D>> passedBox)
      : levelSet(passedLevelSet), box(passedBox) {
    geometry = GeometryEnum::BOX;
  }

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet,
               SmartPointer<Cylinder<T, D>> passedCylinder)
      : levelSet(passedLevelSet), cylinder(passedCylinder) {
    geometry = GeometryEnum::CYLINDER;
  }

  MakeGeometry(SmartPointer<Domain<T, D>> passedLevelSet,
               SmartPointer<PointCloud<T, D>> passedPointCloud)
      : levelSet(passedLevelSet), pointCloud(passedPointCloud) {
    geometry = GeometryEnum::CUSTOM;
  }

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set sphere as geometry to be created in the level set.
  void setGeometry(SmartPointer<Sphere<T, D>> passedSphere) {
    sphere = passedSphere;
    geometry = GeometryEnum::SPHERE;
  }

  /// Set a plane to be created in the level set.
  void setGeometry(SmartPointer<Plane<T, D>> passedPlane) {
    plane = passedPlane;
    geometry = GeometryEnum::PLANE;
  }

  /// Set a box to be created in the level set.
  void setGeometry(SmartPointer<Box<T, D>> passedBox) {
    box = passedBox;
    geometry = GeometryEnum::BOX;
  }

  /// Set a cylinder to be created in the level set.
  void setGeometry(SmartPointer<Cylinder<T, D>> passedCylinder) {
    cylinder = passedCylinder;
    geometry = GeometryEnum::CYLINDER;
  }

  /// Set a point cloud, which is used to create
  /// a geometry from its convex hull.
  void setGeometry(SmartPointer<PointCloud<T, D>> passedPointCloud) {
    pointCloud = passedPointCloud;
    geometry = GeometryEnum::CUSTOM;
  }

  /// Ignore boundary conditions, meaning the parts of the generated
  /// geometry which are outside of the domain boundaries are ignored.
  void setIgnoreBoundaryConditions(bool passedIgnoreBoundaryConditions) {
    for (unsigned i = 0; i < D; ++i) {
      ignoreBoundaryConditions[i] = passedIgnoreBoundaryConditions;
    }
  }

  /// Ignore boundary conditions, meaning the parts of the generated
  /// geometry which are outside of the domain boundaries are ignored.
  /// Set it for each direction separately.
  template <std::size_t N>
  void setIgnoreBoundaryConditions(
      std::array<bool, N> passedIgnoreBoundaryConditions) {
    for (unsigned i = 0; i < D && i < N; ++i) {
      ignoreBoundaryConditions[i] = passedIgnoreBoundaryConditions[i];
    }
  }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to MakeGeometry.")
          .print();
      return;
    }

    switch (geometry) {
    case GeometryEnum::SPHERE:
      makeSphere(sphere->origin, sphere->radius);
      break;
    case GeometryEnum::PLANE:
      makePlane(plane->origin, plane->normal);
      break;
    case GeometryEnum::BOX:
      makeBox(box->minCorner, box->maxCorner);
      break;
    case GeometryEnum::CYLINDER:
      makeCylinder(cylinder);
      break;
    case GeometryEnum::CUSTOM:
      makeCustom(pointCloud);
      break;
    default:
      Logger::getInstance()
          .addWarning("Invalid geometry type was specified for MakeGeometry. "
                      "Not creating geometry.")
          .print();
    }
  }

private:
  void makeSphere(VectorType<T, D> origin, T radius) {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to MakeGeometry.")
          .print();
      return;
    }

    // TODO, this is a stupid algorithm and scales with volume, which is madness
    auto &grid = levelSet->getGrid();
    viennahrle::CoordType gridDelta = grid.getGridDelta();

    // calculate indices from sphere size
    viennahrle::Index<D> index;
    viennahrle::Index<D> endIndex;

    for (unsigned i = 0; i < D; ++i) {
      index[i] = (origin[i] - radius) / gridDelta - 1;
      endIndex[i] = (origin[i] + radius) / gridDelta + 1;
    }

    constexpr double initialWidth = 2.;
    const T valueLimit = initialWidth * 0.5 * gridDelta + 1e-5;
    const T radius2 = radius * radius;

    pointDataType pointData;
    const viennahrle::Index<D> minIndex = index;

    while (index < endIndex) {
      // take the shortest manhattan distance to gridline intersection
      T distance = std::numeric_limits<T>::max();
      for (unsigned i = 0; i < D; ++i) {
        T y = (index[(i + 1) % D] * gridDelta) - origin[(i + 1) % D];
        T z = 0;
        if constexpr (D == 3)
          z = (index[(i + 2) % D] * gridDelta) - origin[(i + 2) % D];
        T x = radius2 - y * y - z * z;
        if (x < 0.)
          continue;
        T dirRadius =
            std::abs((index[i] * gridDelta) - origin[i]) - std::sqrt(x);
        if (std::abs(dirRadius) < std::abs(distance))
          distance = dirRadius;
      }

      if (std::abs(distance) <= valueLimit) {
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
    for (unsigned i = 0; i < pointData.size(); ++i) {
      for (unsigned j = 0; j < D; ++j) {
        if (!ignoreBoundaryConditions[j] && grid.isBoundaryPeriodic(j)) {
          pointData[i].first[j] =
              grid.globalIndex2LocalIndex(j, pointData[i].first[j]);
        }
      }
    }

    levelSet->insertPoints(pointData);
    levelSet->getDomain().segment();
    levelSet->finalize(initialWidth);
  }

  /// Creates a plane containing the point origin, with
  /// the plane normal given by normal
  void makePlane(VectorType<T, D> origin,
                 VectorType<T, D> const &passedNormal) {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to MakeGeometry.")
          .print();
      return;
    }

    auto &grid = levelSet->getGrid();
    viennahrle::CoordType gridDelta = grid.getGridDelta();

    // normalise passedNormal
    double modulus = 0.;
    VectorType<T, D> normal = passedNormal;
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
          viennahrle::BoundaryType::INFINITE_BOUNDARY) {
        if (!infDimSet) {
          i = n;
          infDimSet = true;
        } else {
          Logger::getInstance().addError(
              "Planes can only be created with one Infinite Boundary "
              "Condition. More than one found!");
        }
      }
    }
    if (!infDimSet) {
      Logger::getInstance().addError("Planes require exactly one Infinite "
                                     "Boundary Condition. None found!");
    }

    if (passedNormal[i] == 0.) {
      Logger::getInstance().addError(
          "MakeGeometry: Plane cannot be parallel to Infinite Boundary "
          "direction!");
    }

    // find minimum and maximum points in infinite direction
    // there are 2*(D-1) points in the corners of the simulation domain
    std::vector<Vec3D<T>> cornerPoints;
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
      minCoord[n] = gridDelta * (grid.getMinIndex((i + n + 1) % D) - 1);
      maxCoord[n] = gridDelta * (grid.getMaxIndex((i + n + 1) % D) + 1);
    }

    // set corner points
    cornerPoints[0][j] = minCoord[0];
    cornerPoints[1][j] = maxCoord[0];

    if constexpr (D == 3) {
      cornerPoints[0][k] = minCoord[1];
      cornerPoints[1][k] = maxCoord[1];

      cornerPoints[2][j] = minCoord[0];
      cornerPoints[2][k] = maxCoord[1];
      cornerPoints[3][j] = maxCoord[0];
      cornerPoints[3][k] = minCoord[1];
    }

    // now find i coordinate of points
    auto mesh = SmartPointer<Mesh<T>>::New();

    for (unsigned n = 0; n < cornerPoints.size(); ++n) {
      double numerator = (cornerPoints[n][j] - origin[j]) * normal[j];
      if constexpr (D == 3)
        numerator += (cornerPoints[n][k] - origin[k]) * normal[k];
      else
        cornerPoints[n][2] = 0.;
      cornerPoints[n][i] = origin[i] - numerator / normal[i];
      mesh->insertNextNode(cornerPoints[n]);
    }

    if (D == 2) {
      std::array<unsigned, 2> line = {0, 1};
      if (normal[i] < 0.)
        std::swap(line[0], line[1]);
      mesh->insertNextLine(line);
    } else {
      std::array<unsigned, 3> triangle = {0, 1, 2};
      if (normal[i] < 0.)
        std::swap(triangle[0], triangle[1]);
      mesh->insertNextTriangle(triangle);
      triangle = {0, 3, 1};
      if (normal[i] < 0.)
        std::swap(triangle[0], triangle[1]);
      mesh->insertNextTriangle(triangle);
    }

#ifndef NDEBUG
    static unsigned planeCounter = 0;
    VTKWriter<T>(mesh, "plane" + std::to_string(planeCounter++) + ".vtk")
        .apply();
#endif

    // now convert mesh to levelset
    FromSurfaceMesh<T, D>(levelSet, mesh).apply();
  }

  // This function creates a box starting in minCorner spanning to maxCorner
  void makeBox(VectorType<T, D> minCorner, VectorType<T, D> maxCorner) {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to MakeGeometry.")
          .print();
      return;
    }

    // draw all triangles for the surface and then import from the mesh
    std::vector<Vec3D<T>> corners;
    corners.resize(std::pow(2, D), Vec3D<T>{0, 0, 0});

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

    if constexpr (D == 3) {
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
    auto mesh = Mesh<T>::New();
    for (unsigned i = 0; i < corners.size(); ++i) {
      mesh->insertNextNode(corners[i]);
    }

    if (D == 2) {
      std::array<unsigned, 2> lines[4] = {{0, 2}, {2, 3}, {3, 1}, {1, 0}};
      for (auto &line : lines)
        mesh->insertNextLine(line);
    } else {
      std::array<unsigned, 3> triangles[12] = {
          {0, 3, 1}, {0, 2, 3}, {0, 1, 5}, {0, 5, 4}, {0, 4, 2}, {4, 6, 2},
          {7, 6, 4}, {7, 4, 5}, {7, 2, 6}, {7, 3, 2}, {1, 3, 5}, {3, 7, 5}};
      for (auto &triangle : triangles)
        mesh->insertNextTriangle(triangle);
    }

    // now convert mesh to levelset
    FromSurfaceMesh<T, D> mesher(levelSet, mesh);
    mesher.setRemoveBoundaryTriangles(ignoreBoundaryConditions);
    mesher.apply();
  }

  void makeCylinder(SmartPointer<Cylinder<T, D>> cylinder) {
    if (D != 3) {
      Logger::getInstance()
          .addWarning("MakeGeometry: Cylinder can only be created in 3D!")
          .print();
      return;
    }
    // generate the points on the edges of the cylinders and mesh
    // them manually
    // cylinder axis will be (0,0,1)
    auto gridDelta = levelSet->getGrid().getGridDelta();

    auto points = SmartPointer<PointCloud<T, D>>::New();
    const unsigned numPoints =
        std::ceil(2 * M_PI * cylinder->radius / gridDelta);
    const double smallAngle = 2.0 * M_PI / static_cast<double>(numPoints);

    auto mesh = SmartPointer<Mesh<T>>::New();
    // insert midpoint at base
    mesh->insertNextNode(Vec3D<T>{0.0, 0.0, 0.0});
    {
      constexpr double limit = 2 * M_PI - 1e-6;
      std::vector<Vec3D<T>> points;
      if (cylinder->topRadius)
        std::vector<Vec3D<T>> pointsTop;

      // create and insert points at base
      for (double angle = 0.; angle < limit; angle += smallAngle) {
        Vec3D<T> point;
        point[0] = cylinder->radius * std::cos(angle);
        point[1] = cylinder->radius * std::sin(angle);
        point[2] = 0.0;
        points.push_back(point);
        mesh->insertNextNode(point);
      }

      // insert midpoint at top
      mesh->insertNextNode(Vec3D<T>{0.0, 0.0, cylinder->height});

      double angle = 0;
      for (unsigned i = 0; i < numPoints; ++i) {
        // create triangles at base
        std::array<unsigned, 3> triangle{};
        triangle[0] = (i + 1) % numPoints + 1;
        triangle[1] = i + 1;
        triangle[2] = 0;
        mesh->insertNextTriangle(triangle);

        // insert points at top
        // If topRadius is specified, update the first two coordinates of the
        // points
        if (cylinder->topRadius) {
          points[i][0] = cylinder->topRadius * std::cos(angle);
          points[i][1] = cylinder->topRadius * std::sin(angle);
          angle += smallAngle;
        }
        points[i][2] = cylinder->height;
        mesh->insertNextNode(points[i]);

        // insert triangles at top
        triangle[0] = numPoints + 1;
        triangle[1] = numPoints + i + 2;
        triangle[2] = (i + 1) % numPoints + 2 + numPoints;
        mesh->insertNextTriangle(triangle);
      }

      // insert sidewall triangles
      for (unsigned i = 0; i < numPoints; ++i) {
        std::array<unsigned, 3> triangle{};
        triangle[0] = i + 1;
        triangle[1] = (i + 1) % numPoints + 1;
        triangle[2] = i + numPoints + 2;
        mesh->insertNextTriangle(triangle);

        triangle[0] = (i + 1) % numPoints + 1;
        triangle[1] = (i + 1) % numPoints + 2 + numPoints;
        triangle[2] = i + numPoints + 2;
        mesh->insertNextTriangle(triangle);
      }
    }

    // rotate mesh
    // normalise axis vector
    T unit =
        std::sqrt(DotProduct(cylinder->axisDirection, cylinder->axisDirection));
    Vec3D<T> cylinderAxis;
    for (int i = 0; i < 3; ++i) {
      cylinderAxis[i] = cylinder->axisDirection[i] / unit;
    }
    // get rotation axis via cross product of (0,0,1) and axis of cylinder
    Vec3D<T> rotAxis = {-cylinderAxis[1], cylinderAxis[0], 0.0};
    // angle is acos of dot product
    T rotationAngle = std::acos(cylinderAxis[2]);

    // rotate mesh
    TransformMesh<T>(mesh, TransformEnum::ROTATION, rotAxis, rotationAngle)
        .apply();

    // translate mesh
    Vec3D<T> translationVector;
    for (int i = 0; i < 3; ++i) {
      translationVector[i] = cylinder->origin[i];
    }
    TransformMesh<T>(mesh, TransformEnum::TRANSLATION, translationVector)
        .apply();

    // read mesh from surface
    FromSurfaceMesh<T, D> mesher(levelSet, mesh);
    mesher.setRemoveBoundaryTriangles(ignoreBoundaryConditions);
    mesher.apply();
  }

  void makeCustom(SmartPointer<PointCloud<T, D>> pointCloud) {
    // create mesh from point cloud
    auto mesh = SmartPointer<Mesh<T>>::New();
    ConvexHull<T, D>(mesh, pointCloud).apply();

    // read mesh from surface
    FromSurfaceMesh<T, D> mesher(levelSet, mesh);
    mesher.setRemoveBoundaryTriangles(ignoreBoundaryConditions);
    mesher.apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(MakeGeometry)

} // namespace viennals
