#ifndef LS_BOOLEAN_OPERATION_TEMPLATE_HPP
#define LS_BOOLEAN_OPERATION_TEMPLATE_HPP

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain_template.hpp>

template <class T, int D> class lsBooleanOperation {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;
  lsDomain<T, D> &levelSetA;

  void booleanOpInternal(const lsDomain<T, D> &levelSetB,
                         const T &(*comp)(const T &, const T &)) {
    auto &grid = levelSetA.getGrid();
    lsDomain<T, D> newlsDomain(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain.getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSetA.getDomain();

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

      hrleConstSparseIterator<hrleDomainType> itA(levelSetA.getDomain(),
                                                  currentVector);
      hrleConstSparseIterator<hrleDomainType> itB(levelSetB.getDomain(),
                                                  currentVector);

      while (currentVector < endVector) {

        T currentValue = comp(itA.getValue(), itB.getValue());

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
    levelSetA.deepCopy(newlsDomain);
    levelSetA.finalize(
        std::min(levelSetA.getLevelSetWidth(), levelSetB.getLevelSetWidth()));
  }

  static const T &minComp(const T &A, const T &B) { return std::min(A, B); }

  static const T &maxComp(const T &A, const T &B) { return std::max(A, B); }

  static const T &relativeComplementComp(const T &A, const T &B) {
    return std::max(A, -B);
  }

public:
  lsBooleanOperation(lsDomain<T, D> &passedlsDomain)
      : levelSetA(passedlsDomain){};

  void intersect(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::maxComp);
  }

  void AND(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::maxComp);
  }

  void max(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::maxComp);
  }

  void unite(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::minComp);
  }

  void OR(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::minComp);
  }

  void min(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::minComp);
  }

  void relativeComplement(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::relativeComplementComp);
  }

  void XOR(const lsDomain<T, D> &levelSetB) {
    booleanOpInternal(levelSetB, &lsBooleanOperation::relativeComplementComp);
  }

  void invert() {
    auto &hrleDomain = levelSetA.getDomain();
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
    levelSetA.finalize();
  }

  void NOT() { invert(); }
};

#endif // LS_BOOLEAN_OPERATION_TEMPLATE_HPP
