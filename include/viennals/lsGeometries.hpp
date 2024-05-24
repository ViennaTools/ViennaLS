#pragma once

#include <cassert>
#include <vector>

#include <hrleVectorType.hpp>

#include <lsPreCompileMacros.hpp>

namespace viennals {

/// Class describing a sphere via origin and radius.
template <class T, int D> class Sphere {
public:
  hrleVectorType<T, D> origin = hrleVectorType<T, D>(T(0));
  T radius = 0.;

  Sphere() {}

  Sphere(hrleVectorType<T, D> passedOrigin, T passedRadius)
      : origin(passedOrigin), radius(passedRadius) {}

  Sphere(T *passedOrigin, T passedRadius) : radius(passedRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
    }
  }

  Sphere(const std::vector<T> &passedOrigin, T passedRadius)
      : origin(passedOrigin), radius(passedRadius) {}
};

/// Class describing a plane via a point in it and the plane normal.
template <class T, int D> class Plane {
public:
  hrleVectorType<T, D> origin = hrleVectorType<T, D>(T(0));
  hrleVectorType<T, D> normal = hrleVectorType<T, D>(T(0));

  Plane() {}

  Plane(hrleVectorType<T, D> passedOrigin, hrleVectorType<T, D> passedNormal)
      : origin(passedOrigin), normal(passedNormal) {}

  Plane(const T *passedOrigin, const T *passedNormal) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      normal[i] = passedNormal[i];
    }
  }

  Plane(const std::vector<T> &passedOrigin, const std::vector<T> &passedNormal)
      : origin(passedOrigin), normal(passedNormal) {}
};

/// Class describing a square box from one coordinate to another.
template <class T, int D> class Box {
public:
  hrleVectorType<T, D> minCorner = hrleVectorType<T, D>(T(0));
  hrleVectorType<T, D> maxCorner = hrleVectorType<T, D>(T(0));

  Box() {}

  Box(hrleVectorType<T, D> passedMinCorner,
      hrleVectorType<T, D> passedMaxCorner)
      : minCorner(passedMinCorner), maxCorner(passedMaxCorner) {}

  Box(const T *passedMinCorner, const T *passedMaxCorner) {
    for (unsigned i = 0; i < D; ++i) {
      minCorner[i] = passedMinCorner[i];
      maxCorner[i] = passedMaxCorner[i];
    }
  }

  Box(const std::vector<T> &passedMinCorner,
      const std::vector<T> &passedMaxCorner)
      : minCorner(passedMinCorner), maxCorner(passedMaxCorner) {}
};

/// Class describing a square box from one coordinate to another.
template <class T, int D> class Cylinder {
public:
  /// This is the location of the center of the base of the cylinder
  hrleVectorType<T, 3> origin = hrleVectorType<T, 3>(T(0));
  /// This vector will be the main axis of the cylinder
  hrleVectorType<T, 3> axisDirection = hrleVectorType<T, 3>(T(0));
  /// height of the cylinder
  T height = 0.;
  /// radius of the base of the cylinder
  T radius = 0.;
  /// radius of the top of the cylinder
  T topRadius = 0.;

  Cylinder() {}

  Cylinder(hrleVectorType<T, D> passedOrigin,
           hrleVectorType<T, D> passedAxisDirection, T passedHeight,
           T passedRadius, T passedTopRadius = 0)
      : origin(passedOrigin), axisDirection(passedAxisDirection),
        height(passedHeight), radius(passedRadius), topRadius(passedTopRadius) {
  }

  Cylinder(const T *passedOrigin, const T *passedAxisDirection,
           const T passedHeight, const T passedRadius,
           const T passedTopRadius = 0)
      : height(passedHeight), radius(passedRadius), topRadius(passedTopRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      axisDirection[i] = passedAxisDirection[i];
    }
  }

  Cylinder(std::vector<T> passedOrigin, std::vector<T> passedAxisDirection,
           T passedHeight, T passedRadius, T passedTopRadius = 0)
      : origin(passedOrigin), axisDirection(passedAxisDirection),
        height(passedHeight), radius(passedRadius), topRadius(passedTopRadius) {
  }
};

/// Class describing a point cloud, which can be used to
/// create geometries from its convex hull mesh.
template <class T, int D> class PointCloud {
public:
  std::vector<hrleVectorType<T, D>> points;

  PointCloud() {}

  PointCloud(std::vector<hrleVectorType<T, D>> passedPoints)
      : points(passedPoints) {}

  PointCloud(const std::vector<std::vector<T>> &passedPoints) {
    for (auto point : passedPoints) {
      hrleVectorType<T, D> p(point);
      points.push_back(p);
    }
  }

  void insertNextPoint(hrleVectorType<T, D> newPoint) {
    points.push_back(newPoint);
  }

  void insertNextPoint(T *newPoint) {
    hrleVectorType<T, D> point(newPoint);
    points.push_back(point);
  }

  void insertNextPoint(const std::array<T, D> newPoint) {
    hrleVectorType<T, D> point(newPoint);
    points.push_back(std::move(point));
  }

  void insertNextPoint(const std::vector<T> &newPoint) {
    hrleVectorType<T, D> point(newPoint);
    points.push_back(point);
  }

  std::pair<typename std::vector<hrleVectorType<T, D>>::iterator, bool>
  insertNextUniquePoint(hrleVectorType<T, D> newPoint) {
    for (auto it = points.begin(); it != points.end(); ++it) {
      if (newPoint == *it)
        return std::make_pair(it, false);
    }
    points.push_back(newPoint);
    return std::make_pair(--points.end(), true);
  }

  typename std::vector<hrleVectorType<T, D>>::iterator begin() {
    return points.begin();
  }

  typename std::vector<hrleVectorType<T, D>>::iterator end() {
    return points.end();
  }

  std::size_t size() { return points.size(); }

  hrleVectorType<T, D> &operator[](std::size_t i) { return points[i]; }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Sphere)
PRECOMPILE_PRECISION_DIMENSION(Plane)
PRECOMPILE_PRECISION_DIMENSION(Box)
PRECOMPILE_PRECISION_DIMENSION(PointCloud)

} // namespace viennals
