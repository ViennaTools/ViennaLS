#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>
#include <lsPrune.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

#include <functional>

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

  void booleanOpInternal(ComparatorType comp) {
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

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

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

      // scalars
      for (unsigned i = 0; i < AData.getScalarDataSize(); ++i) {
        auto scalarDataLabel = AData.getScalarDataLabel(i);
        auto BPointer = BData.getScalarData(scalarDataLabel, true);
        if (BPointer != nullptr) {
          auto APointer = AData.getScalarData(i);
          // copy all data into the new scalarData
          typename Domain<T, D>::PointDataType::ScalarDataType scalars;
          scalars.resize(newlsDomain->getNumberOfPoints());
          for (unsigned j = 0; j < newlsDomain->getNumberOfPoints(); ++j) {
            scalars[j] = (newDataLS[0][j])
                             ? APointer->at(newDataSourceIds[0][j])
                             : BPointer->at(newDataSourceIds[0][j]);
          }
          newData.insertNextScalarData(scalars, scalarDataLabel);
        }
      }

      // vectors
      for (unsigned i = 0; i < AData.getVectorDataSize(); ++i) {
        auto vectorDataLabel = AData.getVectorDataLabel(i);
        auto BPointer = BData.getVectorData(vectorDataLabel, true);
        if (BPointer != nullptr) {
          auto APointer = AData.getVectorData(i);
          // copy all data into the new vectorData
          typename Domain<T, D>::PointDataType::VectorDataType vectors;
          vectors.resize(newlsDomain->getNumberOfPoints());
          for (unsigned j = 0; j < newlsDomain->getNumberOfPoints(); ++j) {
            vectors[j] = (newDataLS[0][j])
                             ? APointer->at(newDataSourceIds[0][j])
                             : BPointer->at(newDataSourceIds[0][j]);
          }
          newData.insertNextVectorData(vectors, vectorDataLabel);
        }
      }
    }

    newDomain.finalize();
    newDomain.segment();
    newlsDomain->setLevelSetWidth(levelSetA->getLevelSetWidth());

    if (pruneResult) {
      auto pruner = Prune<T, D>(newlsDomain);
      pruner.setRemoveStrayZeros(true);
      pruner.apply();

      // now we need to prune, to remove stray defined points
      Prune<T, D>(newlsDomain).apply();
    }

    levelSetA->deepCopy(newlsDomain);
  }

  void invert() {
    auto &hrleDomain = levelSetA->getDomain();
#pragma omp parallel num_threads(hrleDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif
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
      booleanOpInternal(&BooleanOperation::maxComp);
      break;
    case BooleanOperationEnum::UNION:
      booleanOpInternal(&BooleanOperation::minComp);
      break;
    case BooleanOperationEnum::RELATIVE_COMPLEMENT:
      booleanOpInternal(&BooleanOperation::relativeComplementComp);
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
