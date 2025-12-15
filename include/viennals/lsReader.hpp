#pragma once

#include <fstream>

#include <lsDomain.hpp>
#include <lsPreCompileMacros.hpp>
#include <utility>

namespace viennals {

using namespace viennacore;

template <class T, int D> class Reader {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  std::string fileName;

public:
  Reader() = default;

  Reader(SmartPointer<Domain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  Reader(SmartPointer<Domain<T, D>> passedLevelSet, std::string passedFileName)
      : levelSet(passedLevelSet), fileName(std::move(passedFileName)) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// set file name for file to write
  void setFileName(std::string passedFileName) {
    fileName = std::move(passedFileName);
  }

  void apply() {
    // check level-set
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level-set was passed to Reader.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      Logger::getInstance()
          .addError("No file name specified for Reader.")
          .print();
      return;
    }
    if (fileName.find(".lvst") != fileName.length() - 5) {
      VIENNACORE_LOG_WARNING(
          "File name does not end in '.lvst', appending it.");
      fileName.append(".lvst");
    }

    // Open file for writing and save serialized level set in it
    std::ifstream fin(fileName, std::ios::binary);

    levelSet->deserialize(fin);

    fin.close();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Reader)

} // namespace viennals
