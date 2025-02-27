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
    // check mesh
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to Reader. Not writing.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      Logger::getInstance()
          .addWarning("No file name specified for Reader. Not writing.")
          .print();
      return;
    }
    if (fileName.find(".lvst") != fileName.length() - 5) {
      Logger::getInstance()
          .addWarning("File name does not end in '.lvst', appending it.")
          .print();
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
