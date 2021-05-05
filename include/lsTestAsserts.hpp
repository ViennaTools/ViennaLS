#pragma once

#include <lsCheck.hpp>

#define LSTEST_ASSERT(condition)                                               \
  {                                                                            \
    if (!(condition)) {                                                        \
      throw std::runtime_error(std::string(__FILE__) + std::string(":") +      \
                               std::to_string(__LINE__) +                      \
                               std::string(" in ") +                           \
                               std::string(__PRETTY_FUNCTION__));              \
    }                                                                          \
  }

#define LSTEST_ASSERT_VALID_LS(levelSet, NumericType, D)                       \
  {                                                                            \
    auto check = lsCheck<NumericType, D>(levelSet);                            \
    check.apply();                                                             \
    if (check.isValid()) {                                                     \
      std::cout << "SUCCESS" << std::endl;                                     \
    } else {                                                                   \
      throw std::runtime_error(                                                \
          std::string(__FILE__) + std::string(":") +                           \
          std::to_string(__LINE__) + std::string(" in ") +                     \
          std::string(__PRETTY_FUNCTION__) + "\n" + check.what());             \
    }                                                                          \
  }
