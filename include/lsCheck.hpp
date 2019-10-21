#ifndef LS_CHECK_HPP
#define LS_CHECK_HPP

#include <ostream>
#include <string>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain_template.hpp>

///  This class is used to find errors in the underlying level set
///  structure, like invalid neighbours of different signs.
template <class T, int D> class lsCheck {
  const lsDomain<T, D> &levelSet;

  int GetStatusFromDistance(T value) {
    int x = static_cast<int>(value);
    if (value >= 0.) {
      return (value <= static_cast<T>(x) + static_cast<T>(0.5)) ? x : x + 1;
    } else {
      return (value >= static_cast<T>(x) - static_cast<T>(0.5)) ? x : x - 1;
    }
  }

public:
  lsCheck(const lsDomain<T, D> &passedLevelSet) : levelSet(passedLevelSet) {}

  std::string apply() {
    std::ostringstream oss;

    for (hrleSparseStarIterator<lsDomain<T, D>> it(levelSet); !it.isFinished();
         it.next()) {

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
                  << it.getCenter().end_indices()
                  << " has undefined neighbor grid points of opposite sign in "
                     "direction "
                  << i << "!" << std::endl;
            }
          }
        }
      }
    }

    return oss.str();
  }
}

#endif // LS_CHECK_HPP
