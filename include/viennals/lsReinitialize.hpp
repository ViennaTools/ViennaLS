#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsToDualContouringMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Exact Euclidean values with their closest points and normals, recomputed
/// from the dual contouring mesh of the planes a level set stores, on a band
/// of the given width. It restores a Euclidean level set after operations
/// that change its values only approximately, such as the minimum or maximum
/// of a Boolean operation near the edges it creates, and widens its band.
/// Sharp edges and corners of the stored planes stay sharp. A level set
/// without stored planes is converted with the planes the dual contouring
/// estimates. Other point data, such as coverages, move to the new band by
/// grid index; the band state of AdvectHermite is dropped.
template <class T, int D> class Reinitialize {
  using hrleIndex = viennahrle::Index<D>;
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  int width = 0;

  // point data that is not Hermite data and moves to the new band
  static bool carried(const std::string &label) {
    return label != FromSurfaceMesh<T, D>::closestPointsLabel &&
           label != FromSurfaceMesh<T, D>::surfaceNormalsLabel &&
           label.rfind("Hermite", 0) != 0;
  }

  // values, sides of the undefined runs and stored normals turned, as
  // BooleanOperation's INVERT
  void invert() {
    auto &domain = levelSet->getDomain();
#pragma omp parallel for
    for (unsigned p = 0; p < domain.getNumberOfSegments(); ++p) {
      auto &segment = domain.getDomainSegment(p);
      for (auto &value : segment.definedValues)
        value = -value;
      // both undefined values present, then the runs swapped between them
      if (segment.undefinedValues.size() < 1)
        segment.undefinedValues.push_back(T(Domain<T, D>::NEG_VALUE));
      if (segment.undefinedValues.size() < 2)
        segment.undefinedValues.push_back(
            (segment.undefinedValues[0] == Domain<T, D>::NEG_VALUE)
                ? T(Domain<T, D>::POS_VALUE)
                : T(Domain<T, D>::NEG_VALUE));
      for (unsigned dim = 0; dim < D; ++dim)
        for (auto &run : segment.runTypes[dim]) {
          if (run == viennahrle::RunTypeValues::UNDEF_PT)
            run = viennahrle::RunTypeValues::UNDEF_PT + 1;
          else if (run == viennahrle::RunTypeValues::UNDEF_PT + 1)
            run = viennahrle::RunTypeValues::UNDEF_PT;
        }
    }
    levelSet->finalize();
    if (auto *normals = levelSet->getPointData().getVectorData(
            FromSurfaceMesh<T, D>::surfaceNormalsLabel, true))
      for (auto &normal : *normals)
        for (auto &component : normal)
          component = -component;
  }

public:
  Reinitialize() = default;

  Reinitialize(SmartPointer<Domain<T, D>> passedLevelSet, int passedWidth = 0)
      : levelSet(passedLevelSet), width(passedWidth) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// Level set width of the new band; 0 keeps the current width.
  void setWidth(int passedWidth) { width = passedWidth; }

  /// Volume (area in 2D) enclosed by a closed mesh, positive where its
  /// normals point out of the material.
  static T orientedVolume(Mesh<T> &mesh) {
    const auto &nodes = mesh.getNodes();
    T volume = 0.;
    for (const auto &e : mesh.template getElements<D>()) {
      if constexpr (D == 2) {
        const auto &a = nodes[e[0]], &b = nodes[e[1]];
        volume += (a[1] * b[0] - a[0] * b[1]) / 2;
      } else {
        const auto &a = nodes[e[0]], &b = nodes[e[1]], &c = nodes[e[2]];
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) -
                   a[1] * (b[0] * c[2] - b[2] * c[0]) +
                   a[2] * (b[0] * c[1] - b[1] * c[0])) /
                  6;
      }
    }
    return volume;
  }

  /// Whether the mesh has no boundary: every node (2D) or edge (3D) used an
  /// even number of times, twice where the mesh is manifold and four times at
  /// the non-manifold edges dual contouring makes where a cell face has
  /// alternating signs.
  static bool isClosed(Mesh<T> &mesh) {
    std::map<std::pair<unsigned, unsigned>, int> uses;
    for (const auto &e : mesh.template getElements<D>()) {
      for (int a = 0; a < D; ++a) {
        if constexpr (D == 2) {
          ++uses[std::make_pair(e[a], e[a])];
        } else {
          const unsigned u = e[a], v = e[(a + 1) % 3];
          ++uses[std::make_pair(std::min(u, v), std::max(u, v))];
        }
      }
    }
    for (const auto &[key, count] : uses)
      if (count % 2 != 0)
        return false;
    return !uses.empty();
  }

  void apply() {
    if (levelSet == nullptr) {
      VIENNACORE_LOG_ERROR("No level set was passed to Reinitialize.");
      return;
    }
    // the mesh of the stored planes with its folds, small flaps facing the
    // wrong way where vertices of neighbouring cells cross over. The sides
    // below take back the nodes a fold turns, and its sharp edges and
    // corners give the distances more accurately than the mesh avoiding
    // folds in 3D
    auto mesh = SmartPointer<Mesh<T>>::New();
    ToDualContouringMesh<T, D> contouring(levelSet, mesh);
    contouring.setAvoidFolds(false);
    contouring.apply();
    // without a surface the values give only the side
    if (mesh->template getElements<D>().empty())
      return;

    // other point data, by grid index, and the old level set for the sides
    auto &pointData = levelSet->getPointData();
    std::map<hrleIndex, unsigned> oldIds;
    for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> it(
             levelSet->getDomain());
         !it.isFinished(); it.next())
      if (it.isDefined())
        oldIds.emplace(it.getStartIndices(), it.getPointId());
    const auto previous = SmartPointer<Domain<T, D>>::New(levelSet);
    std::vector<std::pair<std::string, typename PointData<T>::ScalarDataType>>
        scalars;
    std::vector<std::pair<std::string, typename PointData<T>::VectorDataType>>
        vectors;
    for (unsigned k = 0; k < pointData.getScalarDataSize(); ++k)
      if (carried(pointData.getScalarDataLabel(k)))
        scalars.emplace_back(pointData.getScalarDataLabel(k),
                             *pointData.getScalarData(k));
    for (unsigned k = 0; k < pointData.getVectorDataSize(); ++k)
      if (carried(pointData.getVectorDataLabel(k)))
        vectors.emplace_back(pointData.getVectorDataLabel(k),
                             *pointData.getVectorData(k));

    if (width > 0)
      levelSet->setLevelSetWidth(width);
    levelSet->setNormalization(NormalizationEnum::EUCLIDEAN);
    // the winding number gives the sign for closed surfaces and tolerates
    // non-manifold edges; surfaces ending at the domain boundary take it from
    // the pseudo-normals. A closed surface with the material outside it, as
    // an inverted level set or a void, has negative volume, and the winding
    // number would count nothing as inside: the complement is converted and
    // turned.
    const bool closed = isClosed(*mesh);
    const bool turned = closed && orientedVolume(*mesh) < 0;
    if (turned)
      for (auto &e : mesh->template getElements<D>())
        std::swap(e[0], e[1]);
    FromSurfaceMesh<T, D> mesher(levelSet, mesh);
    mesher.setUseWindingNumber(closed);
    mesher.apply();
    if (turned)
      invert();

    // dual contouring keeps the surface within a cell of where it was, so a
    // node at least one grid spacing from it keeps its side, also a node
    // new to the band, which takes the side of its region; a fold of the
    // mesh, as a small flap facing the wrong way, would turn nodes behind it,
    // since an open surface takes its signs from the pseudo-normals
    {
      auto *storedNormals = levelSet->getPointData().getVectorData(
          FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
      typename Domain<T, D>::PointValueVectorType points;
      bool turnedBack = false;
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> old(
          previous->getDomain());
      for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
               it(levelSet->getDomain());
           !it.isFinished(); it.next()) {
        if (!it.isDefined())
          continue;
        T value = it.getValue();
        old.goToIndicesSequential(it.getStartIndices());
        const T oldValue = old.getValue();
        if (std::abs(oldValue) >= T(1) && value != T(0) &&
            (value < T(0)) != (oldValue < T(0))) {
          value = -value;
          turnedBack = true;
          if (storedNormals != nullptr)
            for (auto &component : (*storedNormals)[it.getPointId()])
              component = -component;
        }
        points.emplace_back(it.getStartIndices(), value);
      }
      // the regions outside the band take their sides from the band again;
      // the same points in the same order keep their point data
      if (turnedBack) {
        const int bandWidth = levelSet->getLevelSetWidth();
        levelSet->insertPoints(points);
        levelSet->getDomain().segment();
        levelSet->finalize(bandWidth);
      }
    }

    // the rebuilt band is a distance field, so neighbours of opposite sign
    // lie within one grid spacing of each other's distance. A node that
    // changed side in the rebuild and breaks this with a neighbour of the
    // other side was turned by a fold of the mesh, where vertices of
    // neighbouring cells cross over, and takes its side back; a turned
    // cluster comes back from its border
    for (int pass = 0; pass < 4; ++pass) {
      auto *storedNormals = levelSet->getPointData().getVectorData(
          FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
      const T slack =
          std::max(T(1e-6), T(1024) * std::numeric_limits<T>::epsilon());
      typename Domain<T, D>::PointValueVectorType points;
      bool turnedBack = false;
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> old(
          previous->getDomain());
      viennahrle::ConstSparseStarIterator<typename Domain<T, D>::DomainType, 1>
          star(levelSet->getDomain());
      for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
               it(levelSet->getDomain());
           !it.isFinished(); it.next()) {
        if (!it.isDefined())
          continue;
        T value = it.getValue();
        old.goToIndicesSequential(it.getStartIndices());
        const T oldValue = old.getValue();
        if (value != T(0) && oldValue != T(0) &&
            (value < T(0)) != (oldValue < T(0))) {
          star.goToIndicesSequential(it.getStartIndices());
          bool broken = false;
          for (int k = 0; k < 2 * D && !broken; ++k) {
            const auto &neighbour = star.getNeighbor(k);
            if (!neighbour.isDefined())
              continue;
            const T other = neighbour.getValue();
            broken = other != T(0) && (other < T(0)) != (value < T(0)) &&
                     std::abs(value) + std::abs(other) > T(1) + slack;
          }
          if (broken) {
            value = -value;
            turnedBack = true;
            if (storedNormals != nullptr)
              for (auto &component : (*storedNormals)[it.getPointId()])
                component = -component;
          }
        }
        points.emplace_back(it.getStartIndices(), value);
      }
      if (!turnedBack)
        break;
      const int bandWidth = levelSet->getLevelSetWidth();
      levelSet->insertPoints(points);
      levelSet->getDomain().segment();
      levelSet->finalize(bandWidth);
    }

    if (scalars.empty() && vectors.empty())
      return;
    // a new node takes the data of the old node at its index; nodes without
    // one take it from their nearest neighbours that have it, layer by layer
    // outwards, as Expand copies point data
    const auto n = levelSet->getNumberOfPoints();
    std::vector<long> source(n, -1);
    std::vector<hrleIndex> indexOf(n);
    std::map<hrleIndex, unsigned> newIds;
    for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> it(
             levelSet->getDomain());
         !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue;
      const auto id = it.getPointId();
      indexOf[id] = it.getStartIndices();
      newIds.emplace(it.getStartIndices(), id);
      auto found = oldIds.find(it.getStartIndices());
      if (found != oldIds.end())
        source[id] = found->second;
    }
    // neighbour offsets, faces first
    std::vector<hrleIndex> offsets;
    for (int c = 0; c < ((D == 2) ? 9 : 27); ++c) {
      hrleIndex o;
      int rest = c, l1 = 0;
      for (int i = 0; i < D; ++i) {
        o[i] = rest % 3 - 1;
        l1 += std::abs(o[i]);
        rest /= 3;
      }
      if (l1 > 0)
        offsets.push_back(o);
    }
    std::stable_sort(offsets.begin(), offsets.end(),
                     [](const hrleIndex &a, const hrleIndex &b) {
                       int la = 0, lb = 0;
                       for (int i = 0; i < D; ++i) {
                         la += std::abs(a[i]);
                         lb += std::abs(b[i]);
                       }
                       return la < lb;
                     });
    for (bool changed = true; changed;) {
      changed = false;
      std::vector<long> next = source;
      for (std::size_t i = 0; i < n; ++i) {
        if (source[i] >= 0)
          continue;
        for (const auto &o : offsets) {
          hrleIndex q = indexOf[i];
          for (int k = 0; k < D; ++k)
            q[k] += o[k];
          auto near = newIds.find(q);
          if (near != newIds.end() && source[near->second] >= 0) {
            next[i] = source[near->second];
            changed = true;
            break;
          }
        }
      }
      source.swap(next);
    }
    for (auto &[label, data] : scalars) {
      typename PointData<T>::ScalarDataType moved(n, T(0));
      for (std::size_t i = 0; i < n; ++i)
        if (source[i] >= 0)
          moved[i] = data[source[i]];
      pointData.insertReplaceScalarData(std::move(moved), label);
    }
    for (auto &[label, data] : vectors) {
      typename PointData<T>::VectorDataType moved(n, Vec3D<T>{});
      for (std::size_t i = 0; i < n; ++i)
        if (source[i] >= 0)
          moved[i] = data[source[i]];
      pointData.insertReplaceVectorData(std::move(moved), label);
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Reinitialize)

} // namespace viennals
