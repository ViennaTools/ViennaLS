#ifndef LS_CHECK_HPP
#define LS_CHECK_HPP

#include <ostream>
#include <string>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>

///  This class is used to find errors in the underlying level set
///  structure, like invalid neighbours of different signs.
template <class T, int D> class lsCheck {
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;

  int GetStatusFromDistance(T value) {
    int x = static_cast<int>(value);
    if (value >= 0.) {
      return (value <= static_cast<T>(x) + static_cast<T>(0.5)) ? x : x + 1;
    } else {
      return (value >= static_cast<T>(x) - static_cast<T>(0.5)) ? x : x - 1;
    }
  }

public:
  lsCheck() {}

  lsCheck(const lsSmartPointer<lsDomain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsCheck.")
          .print();
      return;
    }

    std::ostringstream oss;
    oss << "Report from lsCheck: " << std::endl;

    for (hrleConstSparseStarIterator<hrleDomain<T, D>> it(
             levelSet->getDomain());
         !it.isFinished(); it.next()) {

      if (it.getCenter().isDefined()) {
        for (int i = 0; i < 2 * D; ++i) {
          if (it.getNeighbor(i).isDefined()) {
            if (std::abs(GetStatusFromDistance(it.getCenter().getValue()) -
                         GetStatusFromDistance(it.getNeighbor(i).getValue())) >
                1) {
              oss << "The defined point " << it.getCenter().getStartIndices()
                  << " has an inconsistent defined neighbor in direction " << i
                  << "!" << std::endl;
              oss.precision(24);
              oss << "Value center point: " << it.getCenter().getValue()
                  << "  Value neighbor point: " << it.getNeighbor(i).getValue()
                  << std::endl;
            }
          } else {
            if (it.getNeighbor(i).getValue() >= 0.) { // TODO POS_SIGN
              if (it.getCenter().getValue() < T(-0.5)) {
                oss << "The defined point " << it.getCenter().getStartIndices()
                    << " has a level set value less than -0.5 but has an "
                       "undefined positive neighbor in direction "
                    << i << "!" << std::endl;
              }
            } else {
              if (it.getCenter().getValue() > T(0.5)) {
                oss << "The defined point " << it.getCenter().getStartIndices()
                    << " has a level set value greater than 0.5 but has an "
                       "undefined negative neighbor in direction "
                    << i << "!" << std::endl;
              }
            }
          }
        }
      } else {
        for (int i = 0; i < 2 * D; ++i) {
          if (!it.getNeighbor(i).isDefined()) {
            if ((it.getCenter().getValue() < 0.) !=
                (it.getNeighbor(i).getValue() < 0.)) {
              oss << "The undefined run from "
                  << it.getCenter().getStartIndices() << " to "
                  << it.getCenter().getEndIndices()
                  << " has undefined neighbor grid points of opposite sign in "
                     "direction "
                  << i << "!" << std::endl;
            }
          }
        }
      }
    }

    // output any faults as error
    lsMessage::getInstance().addError(oss.str());
  }
};

#endif // LS_CHECK_HPP
