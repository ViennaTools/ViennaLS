#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsPrune.hpp>
#include <lsReinitialize.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

#include <functional>
#include <memory>
#include <unordered_set>

namespace viennals {

using namespace viennacore;

/// Enumeration for the different types
/// of boolean operations which are
/// supported.
/// When INVERT, only first level set is inverted.
/// When CUSTOM, the user has to supply a valid
/// function pointer of the form const T (*comp)(const T &, const T &).
/// For CUSTOM only the first level set pointer is checked for validity.
enum struct BooleanOperationEnum : unsigned {
  INTERSECT = 0,
  UNION = 1,
  RELATIVE_COMPLEMENT = 2,
  INVERT = 3,
  CUSTOM = 4
};

///  This class is used to perform boolean operations on two
///  level sets and write the resulting level set into the
///  first passed level set.
///  When the boolean operation is set to CUSTOM, a
///  comparator must be set using setBooleanOperationComparator.
///  This comparator returns one value generated from the level set value
///  supplied by each level set. E.g.: for a union, the comparator will
///  always return the smaller of the two values.
///  The function signature for the comparator is defined in the public
///  ComparatorType.
///  When the first level set is Euclidean and an operand stores Hermite data
///  (storesHermiteData), the result is reinitialized
///  (Reinitialize) from the dual contouring mesh of the Hermite data of the
///  operands on the edges the result's surface crosses: on each such edge in
///  each cell the crossing and normal of the operand whose surface bounds the
///  result there, as the operand's own dual contouring mesh takes them in
///  that cell (ToDualContouringMesh::edgeCrossings). Where both operands
///  cross an edge, that is the crossing the result's side changes at: the
///  first one from the end where a maximum is negative or a minimum
///  positive, the last one from the other end, where the value of either
///  end of the edge comes from the second operand beyond the first's
///  (cutTolerance); otherwise the crossing of the first operand, so that an
///  intersection with a superset or a union with a subset leaves the first
///  level set as it is, also where the other's Hermite data come from its
///  values and miss the common surface next to corners. The edges the
///  operation creates stay sharp, also where faces of the operands coincide,
///  as the
///  walls a relative complement removes from a level set they bound
///  together with the other. The minimum and maximum give the sides of the
///  nodes; nodes keep the closest point and normal of the level set their
///  value came from until the reinitialization, and with a CUSTOM comparator
///  the dual contouring mesh of these planes gives the result. INVERT turns
///  the stored normals.
template <class T, int D> class BooleanOperation {
public:
  using ComparatorType =
      std::function<std::pair<T, bool>(const T &, const T &)>;

private:
  typedef typename Domain<T, D>::DomainType hrleDomainType;
  SmartPointer<Domain<T, D>> levelSetA = nullptr;
  SmartPointer<Domain<T, D>> levelSetB = nullptr;
  BooleanOperationEnum operation = BooleanOperationEnum::INTERSECT;
  ComparatorType operationComp = nullptr;
  bool updatePointData = true;
  bool pruneResult = true;

  // how far the value the result takes from the second operand may lie
  // beyond the first's for the first's surface to stay, grid spacings
  static constexpr T cutTolerance = 1e-4;

  struct IndexHash {
    std::size_t operator()(const viennahrle::Index<D> &index) const {
      std::size_t h = 0;
      for (int k = 0; k < D; ++k)
        h = h * 1000003u ^ std::hash<long long>()(index[k]);
      return h;
    }
  };

  // extreme: 1 where the result is the maximum of the operands, -1 where the
  // minimum, 0 for a custom comparator
  void booleanOpInternal(ComparatorType comp, bool negateB = false,
                         int extreme = 0) {
    const bool euclidean =
        levelSetA->getNormalization() == NormalizationEnum::EUCLIDEAN;
    // with Hermite data in an operand; a Euclidean result without them is the
    // minimum or maximum of the values with the pruning of a Manhattan level
    // set (stray zeros where the operands share a face), and keeps its
    // Euclidean normalization (Advect widens its band again)
    const bool hermite = euclidean && (storesHermiteData(levelSetA) ||
                                       storesHermiteData(levelSetB));
    auto &grid = levelSetA->getGrid();
    auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
    typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename Domain<T, D>::DomainType &domain = levelSetA->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    const bool updateData = updatePointData;
    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    std::vector<std::vector<bool>> newDataLS;
    if (updateData) {
      newDataSourceIds.resize(newDomain.getNumberOfSegments());
      newDataLS.resize(newDataSourceIds.size());
    }
    // the nodes whose value comes from B beyond A's (cutTolerance), for the
    // Hermite data of a Euclidean result
    const bool findCuts = hermite && extreme != 0;
    std::vector<std::vector<viennahrle::Index<D>>> cutNodes(
        findCuts ? newDomain.getNumberOfSegments() : 0);

#pragma omp parallel for
    for (unsigned p = 0; p < newDomain.getNumberOfSegments(); ++p) {

      auto &domainSegment = newDomain.getDomainSegment(p);

      viennahrle::Index<D> currentVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      viennahrle::ConstSparseIterator<hrleDomainType> itA(
          levelSetA->getDomain(), currentVector);
      viennahrle::ConstSparseIterator<hrleDomainType> itB(
          levelSetB->getDomain(), currentVector);

      while (currentVector < endVector) {
        const auto &comparison = comp(itA.getValue(), itB.getValue());
        const auto &currentValue = comparison.first;

        if (currentValue != Domain<T, D>::NEG_VALUE &&
            currentValue != Domain<T, D>::POS_VALUE) {
          domainSegment.insertNextDefinedPoint(currentVector, currentValue);
          if (findCuts && !comparison.second &&
              !(std::abs(currentValue - itA.getValue()) <= cutTolerance))
            cutNodes[p].push_back(currentVector);
          if (updateData) {
            // if taken from A, set to true
            const bool originLS = comparison.second;
            newDataLS[p].push_back(originLS);
            const auto originPointId =
                (originLS) ? itA.getPointId() : itB.getPointId();
            newDataSourceIds[p].push_back(originPointId);
          }
        } else {
          domainSegment.insertNextUndefinedPoint(
              currentVector, (currentValue < 0) ? Domain<T, D>::NEG_VALUE
                                                : Domain<T, D>::POS_VALUE);
        }

        switch (Compare(itA.getEndIndices(), itB.getEndIndices())) {
        case -1:
          itA.next();
          break;
        case 0:
          itA.next();
        default:
          itB.next();
        }
        currentVector =
            Max(itA.getStartIndices().get(), itB.getStartIndices().get());
      }
    }

    // merge data
    for (unsigned i = 1; i < newDataLS.size(); ++i) {
      newDataLS[0].insert(newDataLS[0].end(), newDataLS[i].begin(),
                          newDataLS[i].end());

      newDataSourceIds[0].insert(newDataSourceIds[0].end(),
                                 newDataSourceIds[i].begin(),
                                 newDataSourceIds[i].end());
    }

    // transfer data from the old LSs to new LS
    // Only do so if the same data exists in both LSs
    // If this is not the case, the data is invalid
    // and therefore not needed anyway.
    if (updateData) {
      const auto &AData = levelSetA->getPointData();
      const auto &BData = levelSetB->getPointData();
      auto &newData = newlsDomain->getPointData();

      // scalars — preserve A's fields even when B doesn't have them.
      // Previously, fields present only on A (e.g. OxVelocity on the oxide
      // in a RELATIVE_COMPLEMENT with the mask) were silently dropped because
      // the mask has no matching field.  For A-sourced result points we always
      // use A's value; for B-sourced points without a matching B field we use
      // 0.
      for (unsigned i = 0; i < AData.getScalarDataSize(); ++i) {
        auto scalarDataLabel = AData.getScalarDataLabel(i);
        auto APointer = AData.getScalarData(i);
        auto BPointer =
            BData.getScalarData(scalarDataLabel, true); // may be nullptr
        typename Domain<T, D>::PointDataType::ScalarDataType scalars;
        scalars.resize(newlsDomain->getNumberOfPoints(), T(0));
        for (unsigned j = 0; j < newlsDomain->getNumberOfPoints(); ++j) {
          if (newDataLS[0][j])
            scalars[j] = APointer->at(newDataSourceIds[0][j]);
          else if (BPointer != nullptr)
            scalars[j] = BPointer->at(newDataSourceIds[0][j]);
          // else: B-sourced point with no matching B field → stays 0
        }
        newData.insertNextScalarData(scalars, scalarDataLabel);
      }

      // vectors — same asymmetric treatment as scalars above.
      using VectorDataType =
          typename Domain<T, D>::PointDataType::VectorDataType;
      auto transferVectors = [&](const std::string &vectorDataLabel,
                                 const VectorDataType *APointer,
                                 const VectorDataType *BPointer) {
        using Vec = typename VectorDataType::value_type;
        VectorDataType vectors;
        vectors.resize(newlsDomain->getNumberOfPoints(), Vec{});
        // the normals of the second level set turn where its values do
        const bool turn =
            negateB &&
            vectorDataLabel == FromSurfaceMesh<T, D>::surfaceNormalsLabel;
        for (unsigned j = 0; j < newlsDomain->getNumberOfPoints(); ++j) {
          if (newDataLS[0][j]) {
            if (APointer != nullptr)
              vectors[j] = APointer->at(newDataSourceIds[0][j]);
          } else if (BPointer != nullptr) {
            vectors[j] = BPointer->at(newDataSourceIds[0][j]);
            if (turn)
              for (auto &component : vectors[j])
                component = -component;
          }
          // else: B-sourced point with no matching B field → stays
          // zero-initialised
        }
        newData.insertNextVectorData(vectors, vectorDataLabel);
      };
      for (unsigned i = 0; i < AData.getVectorDataSize(); ++i) {
        auto vectorDataLabel = AData.getVectorDataLabel(i);
        transferVectors(vectorDataLabel, AData.getVectorData(i),
                        BData.getVectorData(vectorDataLabel, true));
      }
      // a Euclidean result needs the closest points and normals of B also
      // where A has none, as when A is still empty
      if (hermite)
        for (const auto &label : {FromSurfaceMesh<T, D>::closestPointsLabel,
                                  FromSurfaceMesh<T, D>::surfaceNormalsLabel})
          if (AData.getVectorData(label, true) == nullptr &&
              BData.getVectorData(label, true) != nullptr)
            transferVectors(label, nullptr, BData.getVectorData(label, true));
    }

    newDomain.finalize();
    newDomain.segment();
    newlsDomain->setLevelSetWidth(levelSetA->getLevelSetWidth());
    newlsDomain->setNormalization(levelSetA->getNormalization());

    // a Euclidean result with Hermite data gets its whole band again from its
    // planes
    if (pruneResult && !hermite) {
      auto pruner = Prune<T, D>(newlsDomain);
      pruner.setRemoveStrayZeros(true);
      pruner.apply();

      // now we need to prune, to remove stray defined points
      Prune<T, D>(newlsDomain).apply();
    }

    // the Hermite data of the edges the operands cross, before A changes
    typename ToDualContouringMesh<T, D>::EdgeData edgeData;
    if (hermite && extreme != 0) {
      using Contouring = ToDualContouringMesh<T, D>;
      const auto crossingsA =
          std::make_shared<typename Contouring::EdgeCrossings>(
              Contouring(levelSetA, nullptr).edgeCrossings());
      const auto crossingsB =
          std::make_shared<typename Contouring::EdgeCrossings>(
              Contouring(levelSetB, nullptr).edgeCrossings());
      const auto cuts = std::make_shared<
          std::unordered_set<viennahrle::Index<D>, IndexHash>>();
      for (const auto &nodes : cutNodes)
        cuts->insert(nodes.begin(), nodes.end());
      edgeData = [crossingsA, crossingsB, cuts, negateB, extreme](
                     const viennahrle::Index<D> &cell, int a, int b, T valueA,
                     T, typename Contouring::EdgeCrossing &crossing) {
        int d = 0;
        while (d + 1 < D && (1 << d) != b - a)
          ++d;
        const typename Contouring::EdgeKey key(cell, a * D + d);
        const auto foundA = crossingsA->find(key);
        const auto foundB = crossingsB->find(key);
        const bool hasA = foundA != crossingsA->end(),
                   hasB = foundB != crossingsB->end();
        if (!hasA && !hasB)
          return false;
        bool fromA = hasA;
        // B bounds the result on this edge only where it cuts into A at an
        // end of the edge
        bool cut = false;
        for (const int corner : {a, b}) {
          viennahrle::Index<D> node = cell;
          for (int k = 0; k < D; ++k)
            node[k] += (corner >> k) & 1;
          cut = cut || cuts->count(node) != 0;
        }
        if (hasA && hasB && cut) {
          // along the edge from corner a, which lies at the cell's lowest
          // index along d: a maximum negative at a turns positive at the
          // first crossing, a minimum positive at a negative at the first
          // one, and otherwise the result changes side at the last one
          const T tA = foundA->second.point[d] - T(cell[d]),
                  tB = foundB->second.point[d] - T(cell[d]);
          const bool first = (valueA < T(0)) == (extreme > 0);
          fromA = first ? tA <= tB : tA >= tB;
        }
        crossing = fromA ? foundA->second : foundB->second;
        if (!fromA && negateB)
          for (auto &component : crossing.normal)
            component = -component;
        return true;
      };
    }

    levelSetA->deepCopy(newlsDomain);

    if (hermite) {
      // the band state and plane priorities of AdvectHermite no longer
      // describe the result
      auto &data = levelSetA->getPointData();
      for (int k = int(data.getScalarDataSize()) - 1; k >= 0; --k)
        if (data.getScalarDataLabel(k).rfind("Hermite", 0) == 0)
          data.eraseScalarData(k);
      for (int k = int(data.getVectorDataSize()) - 1; k >= 0; --k)
        if (data.getVectorDataLabel(k).rfind("Hermite", 0) == 0)
          data.eraseVectorData(k);
      Reinitialize<T, D> reinitialize(levelSetA);
      reinitialize.setEdgeData(edgeData);
      reinitialize.apply();
    }
  }

  void invert() {
    auto &hrleDomain = levelSetA->getDomain();
#pragma omp parallel for
    for (unsigned p = 0; p < hrleDomain.getNumberOfSegments(); ++p) {
      auto &domainSegment = hrleDomain.getDomainSegment(p);

      // change all defined values
      for (unsigned i = 0; i < domainSegment.definedValues.size(); ++i) {
        domainSegment.definedValues[i] = -domainSegment.definedValues[i];
      }

      // add undefined values if missing
      if (domainSegment.undefinedValues.size() < 1) {
        domainSegment.undefinedValues.push_back(T(Domain<T, D>::NEG_VALUE));
      }
      if (domainSegment.undefinedValues.size() < 2) {
        if (domainSegment.undefinedValues[0] == Domain<T, D>::NEG_VALUE) {
          domainSegment.undefinedValues.push_back(T(Domain<T, D>::POS_VALUE));
        } else {
          domainSegment.undefinedValues.push_back(T(Domain<T, D>::NEG_VALUE));
        }
      }

      // change all runTypes
      // there are only two undefined runs: negative undefined (UNDEF_PT)
      // and positive undefined (UNDEF_PT+1)
      for (unsigned dim = 0; dim < D; ++dim) {
        for (unsigned c = 0; c < domainSegment.runTypes[dim].size(); ++c) {
          if (domainSegment.runTypes[dim][c] ==
              viennahrle::RunTypeValues::UNDEF_PT) {
            domainSegment.runTypes[dim][c] =
                viennahrle::RunTypeValues::UNDEF_PT + 1;
          } else if (domainSegment.runTypes[dim][c] ==
                     viennahrle::RunTypeValues::UNDEF_PT + 1) {
            domainSegment.runTypes[dim][c] =
                viennahrle::RunTypeValues::UNDEF_PT;
          }
        }
      }
    }
    levelSetA->finalize();

    // closest points stay, the stored normals turn with the values
    if (auto *normals = levelSetA->getPointData().getVectorData(
            FromSurfaceMesh<T, D>::surfaceNormalsLabel, true))
      for (auto &normal : *normals)
        for (auto &component : normal)
          component = -component;
  }

  static std::pair<T, bool> minComp(const T &a, const T &b) {
    bool AIsSmaller = a < b;
    if (AIsSmaller)
      return std::make_pair(a, true);
    else
      return std::make_pair(b, false);
  }

  static std::pair<T, bool> maxComp(const T &a, const T &b) {
    bool AIsLarger = a > b;
    if (AIsLarger)
      return std::make_pair(a, true);
    else
      return std::make_pair(b, false);
  }

  static std::pair<T, bool> relativeComplementComp(const T &a, const T &b) {
    return maxComp(a, -b);
  }

public:
  BooleanOperation() = default;

  BooleanOperation(
      SmartPointer<Domain<T, D>> passedlsDomain,
      BooleanOperationEnum passedOperation = BooleanOperationEnum::INVERT)
      : levelSetA(passedlsDomain), operation(passedOperation) {}

  BooleanOperation(
      SmartPointer<Domain<T, D>> passedlsDomainA,
      SmartPointer<Domain<T, D>> passedlsDomainB,
      BooleanOperationEnum passedOperation = BooleanOperationEnum::INTERSECT)
      : levelSetA(passedlsDomainA), levelSetB(passedlsDomainB),
        operation(passedOperation){};

  /// Set which level set to perform the boolean operation on.
  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSetA = passedlsDomain;
  }

  /// Set the level set which will be used to modify the
  /// first level set.
  void setSecondLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSetB = passedlsDomain;
  }

  /// Set which of the operations of BooleanOperationEnum to perform.
  void setBooleanOperation(BooleanOperationEnum passedOperation) {
    operation = passedOperation;
  }

  /// Set the comparator to be used when the BooleanOperation
  /// is set to CUSTOM.
  void setBooleanOperationComparator(ComparatorType passedOperationComp) {
    operationComp = passedOperationComp;
  }

  /// Set whether to update the point data stored in the LS
  /// during this algorithm. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  /// Set whether the resulting level set should be pruned. Defaults to true
  void setPruneResult(bool pR) { pruneResult = pR; }

  /// Perform operation.
  void apply() {
    if (levelSetA == nullptr) {
      VIENNACORE_LOG_ERROR("No level set was passed to BooleanOperation.");
      return;
    }

    if (static_cast<unsigned>(operation) < 3) {
      if (levelSetB == nullptr) {
        VIENNACORE_LOG_ERROR(
            "Only one level set was passed to BooleanOperation, "
            "although two were required.");
        return;
      }
    }

    switch (operation) {
    case BooleanOperationEnum::INTERSECT:
      booleanOpInternal(&BooleanOperation::maxComp, false, 1);
      break;
    case BooleanOperationEnum::UNION:
      booleanOpInternal(&BooleanOperation::minComp, false, -1);
      break;
    case BooleanOperationEnum::RELATIVE_COMPLEMENT:
      booleanOpInternal(&BooleanOperation::relativeComplementComp, true, 1);
      break;
    case BooleanOperationEnum::INVERT:
      invert();
      break;
    case BooleanOperationEnum::CUSTOM:
      if (operationComp == nullptr) {
        VIENNACORE_LOG_ERROR(
            "No comparator supplied to custom BooleanOperation.");
        return;
      }
      booleanOpInternal(operationComp);
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(BooleanOperation)

} // namespace viennals
