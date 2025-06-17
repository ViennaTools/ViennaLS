#pragma once

#include <cassert>
#include <vector>

#include <lsPreCompileMacros.hpp>

#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {
using namespace viennacore;

/// Class describing a sphere via origin and radius.
template <class T, int D> class Sphere {
public:
  VectorType<T, D> origin;
  T radius = 0.;

  Sphere(VectorType<T, D> passedOrigin, T passedRadius)
      : origin(passedOrigin), radius(passedRadius) {}

  Sphere(T *passedOrigin, T passedRadius) : radius(passedRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
    }
  }

  Sphere(const std::vector<T> &passedOrigin, T passedRadius)
      : radius(passedRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
    }
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<Sphere>::New(std::forward<Args>(args)...);
  }
};

/// Class describing a plane via a point in it and the plane normal.
template <class T, int D> class Plane {
public:
  VectorType<T, D> origin;
  VectorType<T, D> normal;

  Plane(VectorType<T, D> passedOrigin, VectorType<T, D> passedNormal)
      : origin(passedOrigin), normal(passedNormal) {}

  Plane(const T *passedOrigin, const T *passedNormal) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      normal[i] = passedNormal[i];
    }
  }

  Plane(const std::vector<T> &passedOrigin,
        const std::vector<T> &passedNormal) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      normal[i] = passedNormal[i];
    }
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<Plane>::New(std::forward<Args>(args)...);
  }
};

/// Class describing a square box from one coordinate to another.
template <class T, int D> class Box {
public:
  VectorType<T, D> minCorner;
  VectorType<T, D> maxCorner;

  Box(VectorType<T, D> passedMinCorner, VectorType<T, D> passedMaxCorner)
      : minCorner(passedMinCorner), maxCorner(passedMaxCorner) {}

  Box(const T *passedMinCorner, const T *passedMaxCorner) {
    for (unsigned i = 0; i < D; ++i) {
      minCorner[i] = passedMinCorner[i];
      maxCorner[i] = passedMaxCorner[i];
    }
  }

  Box(const std::vector<T> &passedMinCorner,
      const std::vector<T> &passedMaxCorner) {
    for (unsigned i = 0; i < D; ++i) {
      minCorner[i] = passedMinCorner[i];
      maxCorner[i] = passedMaxCorner[i];
    }
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<Box>::New(std::forward<Args>(args)...);
  }
};

/// Class describing a square box from one coordinate to another.
template <class T, int D> class Cylinder {
public:
  /// This is the location of the center of the base of the cylinder
  VectorType<T, 3> origin;
  /// This vector will be the main axis of the cylinder
  VectorType<T, 3> axisDirection;
  /// height of the cylinder
  T height = 0.;
  /// radius of the base of the cylinder
  T radius = 0.;
  /// radius of the top of the cylinder
  T topRadius = 0.;

  Cylinder(VectorType<T, D> passedOrigin, VectorType<T, D> passedAxisDirection,
           T passedHeight, T passedRadius, T passedTopRadius = 0)
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
      : height(passedHeight), radius(passedRadius), topRadius(passedTopRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      axisDirection[i] = passedAxisDirection[i];
    }
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<Cylinder>::New(std::forward<Args>(args)...);
  }
};

/// Class describing a point cloud, which can be used to
/// create geometries from its convex hull mesh.
template <class T, int D> class PointCloud {
public:
  std::vector<VectorType<T, D>> points;

  PointCloud() = default;

  PointCloud(std::vector<VectorType<T, D>> passedPoints)
      : points(passedPoints) {}

  void insertNextPoint(T *newPoint) {
    VectorType<T, D> point;
    for (unsigned i = 0; i < D; ++i) {
      point[i] = newPoint[i];
    }
    points.push_back(std::move(point));
  }

  void insertNextPoint(const VectorType<T, D> &newPoint) {
    VectorType<T, D> point(newPoint);
    points.push_back(std::move(point));
  }

  std::pair<typename std::vector<VectorType<T, D>>::iterator, bool>
  insertNextUniquePoint(VectorType<T, D> newPoint) {
    for (auto it = points.begin(); it != points.end(); ++it) {
      if (newPoint == *it)
        return std::make_pair(it, false);
    }
    points.push_back(newPoint);
    return std::make_pair(--points.end(), true);
  }

  typename std::vector<VectorType<T, D>>::iterator begin() {
    return points.begin();
  }

  typename std::vector<VectorType<T, D>>::iterator end() {
    return points.end();
  }

  std::size_t size() { return points.size(); }

  VectorType<T, D> &operator[](std::size_t i) { return points[i]; }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<PointCloud>::New(std::forward<Args>(args)...);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Sphere)
PRECOMPILE_PRECISION_DIMENSION(Plane)
PRECOMPILE_PRECISION_DIMENSION(Box)
PRECOMPILE_PRECISION_DIMENSION(PointCloud)

} // namespace viennals
