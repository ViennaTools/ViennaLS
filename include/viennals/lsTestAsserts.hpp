#pragma once

#include <cmath>
#include <iostream>
#include <lsCheck.hpp>
#include <lsMesh.hpp>
#include <sstream>
#include <vcTestAsserts.hpp>
#include <vector>

#define LSTEST_ASSERT_VALID_LS(levelSet, NumericType, D)                       \
  {                                                                            \
    auto check = viennals::Check<NumericType, D>(levelSet);                    \
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

#define LSTEST_ASSERT_MESH_CORNERS(mesh, expected, D)                          \
  {                                                                            \
    try {                                                                      \
      auto &nodes = mesh->getNodes();                                          \
      double tolerance = 1e-4;                                                 \
      for (const auto &ex : expected) {                                        \
        bool found = false;                                                    \
        for (const auto &n : nodes) {                                          \
          double distSq = 0.;                                                  \
          for (int i = 0; i < D; ++i)                                          \
            distSq += (n[i] - ex[i]) * (n[i] - ex[i]);                         \
          if (distSq < tolerance) {                                            \
            found = true;                                                      \
            break;                                                             \
          }                                                                    \
        }                                                                      \
        if (!found) {                                                          \
          std::stringstream ss;                                                \
          ss << "Corner not found: (";                                         \
          for (int i = 0; i < D; ++i)                                          \
            ss << ex[i] << (i < D - 1 ? ", " : "");                            \
          ss << ")";                                                           \
          throw std::runtime_error(ss.str());                                  \
        }                                                                      \
      }                                                                        \
      std::cout << "All expected corners found!" << std::endl;                 \
    } catch (const std::exception &e) {                                        \
      throw std::runtime_error(                                                \
          std::string(__FILE__) + std::string(":") +                           \
          std::to_string(__LINE__) + std::string(" in ") +                     \
          std::string(__PRETTY_FUNCTION__) + "\n" + e.what());                 \
    }                                                                          \
  }
