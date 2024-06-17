#pragma once

#include <ostream>
#include <string>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>

namespace viennals {

using namespace viennacore;

enum struct CheckStatusEnum : unsigned {
  SUCCESS = 0,
  FAILED = 1,
  UNCHECKED = 2
};

///  This class is used to find errors in the underlying level set
///  structure, like invalid neighbours of different signs.
template <class T, int D> class Check {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  CheckStatusEnum status = CheckStatusEnum::UNCHECKED;
  std::string errors = "Level Set has not been checked yet!";
  bool printMessage = false;

  std::string getDirectionString(int i) {
    std::string result = "";
    if (i < D) {
      result += "-";
    } else {
      i -= D;
      result += "+";
    }
    if (i == 0) {
      result += "x";
    } else if (i == 1) {
      result += "y";
    } else if (i == 2) {
      result += "z";
    } else {
      Logger::getInstance().addError("Check: " + std::to_string(i) +
                                     " is not a valid direction index!");
    }
    return result;
  }

  bool isNegative(T val) { return val <= -std::numeric_limits<T>::epsilon(); }

  int GetStatusFromDistance(T value) {
    int x = static_cast<int>(value);
    if (value >= 0.) {
      return (value <= static_cast<T>(x) + static_cast<T>(0.5)) ? x : x + 1;
    } else {
      return (value >= static_cast<T>(x) - static_cast<T>(0.5)) ? x : x - 1;
    }
  }

public:
  Check() {}

  Check(SmartPointer<Domain<T, D>> passedLevelSet, bool print = false)
      : levelSet(passedLevelSet), printMessage(print) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setPrintMessage(bool print) { printMessage = print; }

  CheckStatusEnum getStatus() const { return status; }

  bool isValid() const { return status == CheckStatusEnum::SUCCESS; }

  std::string what() const { return errors; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to Check.")
          .print();
      return;
    }

    std::ostringstream oss;

    for (hrleConstSparseStarIterator<hrleDomain<T, D>, 1> it(
             levelSet->getDomain());
         !it.isFinished(); it.next()) {

      if (it.getCenter().isDefined()) {
        for (int i = 0; i < 2 * D; ++i) {
          if (it.getNeighbor(i).isDefined()) {
            if (std::abs(GetStatusFromDistance(it.getCenter().getValue()) -
                         GetStatusFromDistance(it.getNeighbor(i).getValue())) >
                1) {
              oss << "The defined point " << it.getCenter().getStartIndices()
                  << " has an inconsistent defined neighbor in direction "
                  << getDirectionString(i) << "!" << std::endl;
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
                    << getDirectionString(i) << "!" << std::endl;
              }
            } else {
              if (it.getCenter().getValue() > T(0.5)) {
                oss << "The defined point " << it.getCenter().getStartIndices()
                    << " has a level set value greater than 0.5 but has an "
                       "undefined negative neighbor in direction "
                    << getDirectionString(i) << "!" << std::endl;
              }
            }
          }
        }
      } else {
        for (int i = 0; i < 2 * D; ++i) {
          if (!it.getNeighbor(i).isDefined()) {
            if (isNegative(it.getCenter().getValue()) !=
                isNegative(it.getNeighbor(i).getValue())) {
              oss << "The undefined run from "
                  << it.getCenter().getStartIndices() << " to "
                  << it.getCenter().getEndIndices()
                  << " has undefined neighbor grid points of opposite sign in "
                     "direction "
                  << getDirectionString(i) << "!" << std::endl;
            }
          }
        }
      }
    }

    // output any faults as error
    if (std::string s = oss.str(); !s.empty()) {
      status = CheckStatusEnum::FAILED;
      errors = s;
      if (printMessage) {
        std::string message = "Report from Check:\n" + s;
        Logger::getInstance().addError(s);
      }
    } else {
      status = CheckStatusEnum::SUCCESS;
      errors = "";
    }
  }
};

} // namespace viennals
