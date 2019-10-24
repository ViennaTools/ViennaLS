#ifndef LS_BOOLEAN_OPERATION_HPP
#define LS_BOOLEAN_OPERATION_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain.hpp>
#include <lsMessage.hpp>

/// Enumeration for the different types
/// of boolean operations which are
/// supported.
/// WHEN INVERT, only first level set is inverted
/// When CUSTOM, the user has to supply a valid
/// function pointer of the form const T (*comp)(const T &, const T &).
/// For CUSTOM only the first level set pointer is checked for validity.
enum struct lsBooleanOperationEnum : unsigned {
  INTERSECT = 0,
  UNION = 1,
  RELATIVE_COMPLEMENT = 2,
  INVERT = 3,
  CUSTOM = 4
};

///  This class is used to perform boolean operations on two
///  level sets and write the resulting level set into the
///  first passed level set.
template <class T, int D> class lsBooleanOperation {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;
  lsDomain<T, D> *levelSetA = nullptr;
  lsDomain<T, D> *levelSetB = nullptr;
  lsBooleanOperationEnum operation = lsBooleanOperationEnum::INTERSECT;
  const T (*operationComp)(const T &, const T &) = nullptr;

  void booleanOpInternal(const T (*comp)(const T &, const T &)) {
    auto &grid = levelSetA->getGrid();
    lsDomain<T, D> newlsDomain(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain.getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSetA->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &domainSegment = newDomain.getDomainSegment(p);

      hrleVectorType<hrleIndexType, D> currentVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      hrleConstSparseIterator<hrleDomainType> itA(levelSetA->getDomain(),
                                                  currentVector);
      hrleConstSparseIterator<hrleDomainType> itB(levelSetB->getDomain(),
                                                  currentVector);

      while (currentVector < endVector) {

        const T &currentValue = comp(itA.getValue(), itB.getValue());

        if (currentValue != lsDomain<T, D>::NEG_VALUE &&
            currentValue != lsDomain<T, D>::POS_VALUE) {
          domainSegment.insertNextDefinedPoint(currentVector, currentValue);
        } else {
          domainSegment.insertNextUndefinedPoint(
              currentVector, (currentValue < 0) ? lsDomain<T, D>::NEG_VALUE
                                                : lsDomain<T, D>::POS_VALUE);
        }

        switch (compare(itA.getEndIndices(), itB.getEndIndices())) {
        case -1:
          itA.next();
          break;
        case 0:
          itA.next();
        default:
          itB.next();
        }
        currentVector = std::max(itA.getStartIndices(), itB.getStartIndices());
      }
    }
    newDomain.finalize();
    newDomain.segment();
    levelSetA->deepCopy(newlsDomain);
    levelSetA->finalize(
        std::min(levelSetA->getLevelSetWidth(), levelSetB->getLevelSetWidth()));
  }

  static const T minComp(const T &A, const T &B) { return std::min(A, B); }

  static const T maxComp(const T &A, const T &B) { return std::max(A, B); }

  static const T relativeComplementComp(const T &A, const T &B) {
    return std::max(A, -B);
  }

public:
  lsBooleanOperation(lsDomain<T, D> &passedlsDomain,
                     lsBooleanOperationEnum passedOperation =
                         lsBooleanOperationEnum::INTERSECT)
      : levelSetA(&passedlsDomain), operation(passedOperation){};

  lsBooleanOperation(lsDomain<T, D> &passedlsDomainA,
                     lsDomain<T, D> &passedlsDomainB,
                     lsBooleanOperationEnum passedOperation =
                         lsBooleanOperationEnum::INTERSECT)
      : levelSetA(&passedlsDomainA), levelSetB(&passedlsDomainB),
        operation(passedOperation){};

  void setLevelSet(lsDomain<T, D> &passedlsDomain) {
    levelSetA = &passedlsDomain;
  }

  void setSecondLevelSet(lsDomain<T, D> &passedlsDomain) {
    levelSetB = &passedlsDomain;
  }

  void setBooleanOperation(lsBooleanOperationEnum passedOperation) {
    operation = passedOperation;
  }

  void setBooleanOperationComparator(const T (*passedOperationComp)(const T &, const T &)) {
    operationComp = passedOperationComp;
  }

  void apply() {
    if(levelSetA == nullptr) {
      lsMessage::getInstance().addWarning("No level set was passed to lsBooleanOperation. Not performing operation.").print();
      return;
    }

    if(static_cast<unsigned>(operation) < 3) {
      if(levelSetB == nullptr) {
        lsMessage::getInstance().addWarning("Only one level set was passed to lsBooleanOperation, although two were required. Not performing operation.").print();
        return;
      }
    }

    switch(operation){
      case lsBooleanOperationEnum::INTERSECT:
        booleanOpInternal(&lsBooleanOperation::maxComp);
        break;
      case lsBooleanOperationEnum::UNION:
        booleanOpInternal(&lsBooleanOperation::minComp);
        break;
      case lsBooleanOperationEnum::RELATIVE_COMPLEMENT:
        booleanOpInternal(&lsBooleanOperation::relativeComplementComp);
        break;
      case lsBooleanOperationEnum::INVERT:
        invert();
        break;
      case lsBooleanOperationEnum::CUSTOM:
        if(operationComp == nullptr) {
          lsMessage::getInstance().addWarning("No comparator supplied to custom lsBooleanOperation. Not performing operation.").print();
        }
        booleanOpInternal(operationComp);
    }
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
        domainSegment.undefinedValues.push_back(T(lsDomain<T, D>::NEG_VALUE));
      }
      if (domainSegment.undefinedValues.size() < 2) {
        if (domainSegment.undefinedValues[0] == lsDomain<T, D>::NEG_VALUE) {
          domainSegment.undefinedValues.push_back(T(lsDomain<T, D>::POS_VALUE));
        } else {
          domainSegment.undefinedValues.push_back(T(lsDomain<T, D>::NEG_VALUE));
        }
      }

      // change all runTypes
      // there are only two undefined runs: negative undefined (UNDEF_PT)
      // and positive undefined (UNDEF_PT+1)
      for (unsigned dim = 0; dim < D; ++dim) {
        for (unsigned c = 0; c < domainSegment.runTypes[dim].size(); c++) {
          if (domainSegment.runTypes[dim][c] == hrleRunTypeValues::UNDEF_PT) {
            domainSegment.runTypes[dim][c] = hrleRunTypeValues::UNDEF_PT + 1;
          } else if (domainSegment.runTypes[dim][c] ==
                     hrleRunTypeValues::UNDEF_PT + 1) {
            domainSegment.runTypes[dim][c] = hrleRunTypeValues::UNDEF_PT;
          }
        }
      }
    }
    levelSetA->finalize();
  }

  void NOT() { invert(); }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsBooleanOperation)

#endif // LS_BOOLEAN_OPERATION_HPP
