#ifndef LS_GEOMETRIES_HPP
#define LS_GEOMETRIES_HPP

/// Class describing a sphere via origin and radius.
template <class T, int D> class lsSphere {
public:
  hrleVectorType<T, D> origin = hrleVectorType<T, D>(T(0));
  T radius = 0.;

  lsSphere(hrleVectorType<T, D> passedOrigin, T passedRadius)
      : origin(passedOrigin), radius(passedRadius) {}

  lsSphere(T *passedOrigin, T passedRadius) : radius(passedRadius) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
    }
  }

  lsSphere(T x, T y, T z, T passedRadius) : origin(hrleVectorType<T, D>(x, y)), radius(passedRadius) {
    if(D==3) origin[2] = z;
  }

  lsSphere(T x, T y, T passedRadius) : origin(hrleVectorType<T, D>(x, y)), radius(passedRadius) {}
};

/// Class describing a plane via a point in it and the plane normal.
template <class T, int D> class lsPlane {
public:
  hrleVectorType<T, D> origin = hrleVectorType<T, D>(T(0));
  hrleVectorType<T, D> normal = hrleVectorType<T, D>(T(0));

  lsPlane(hrleVectorType<T, D> passedOrigin, hrleVectorType<T, D> passedNormal)
      : origin(passedOrigin), normal(passedNormal) {}

  lsPlane(T *passedOrigin, T *passedNormal) {
    for (unsigned i = 0; i < D; ++i) {
      origin[i] = passedOrigin[i];
      normal[i] = passedNormal[i];
    }
  }
};

/// Class describing a square box from one coordinate to another.
template <class T, int D> class lsBox {
public:
  hrleVectorType<T, D> minCorner = hrleVectorType<T, D>(T(0));
  hrleVectorType<T, D> maxCorner = hrleVectorType<T, D>(T(0));

  lsBox(hrleVectorType<T, D> passedMinCorner,
        hrleVectorType<T, D> passedMaxCorner)
      : minCorner(passedMinCorner), maxCorner(passedMaxCorner) {}

  lsBox(T *passedMinCorner, T *passedMaxCorner) {
    for (unsigned i = 0; i < D; ++i) {
      minCorner[i] = passedMinCorner[i];
      maxCorner[i] = passedMaxCorner[i];
    }
  }
};

/// Class describing a point cloud, which can be used to
/// create geometries from its convex hull mesh.
template <class T, int D> class lsPointCloud {
public:
  std::vector<hrleVectorType<T, D>> points;

  lsPointCloud() {}

  lsPointCloud(std::vector<hrleVectorType<T, D>> passedPoints)
      : points(passedPoints) {}

  void insertNextPoint(hrleVectorType<T, D> newPoint) {
    points.push_back(newPoint);
  }

  void insertNextPoint(T *newPoint) {
    hrleVectorType<T, D> point(newPoint);
    points.push_back(point);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsSphere)
PRECOMPILE_PRECISION_DIMENSION(lsPlane)
PRECOMPILE_PRECISION_DIMENSION(lsBox)
PRECOMPILE_PRECISION_DIMENSION(lsPointCloud)

#endif // LS_GEOMETRIES_HPP
