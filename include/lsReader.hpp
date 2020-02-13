#ifndef LS_READER_HPP
#define LS_READER_HPP

#include <fstream>

#include <lsPreCompileMacros.hpp>
#include <lsDomain.hpp>

template <class T, int D> class lsReader {
  lsDomain<T, D> *levelSet = nullptr;
  std::string fileName;

public:
  lsReader() {}

  lsReader(lsDomain<T, D> &passedLevelSet) : levelSet(&passedLevelSet) {}

  lsReader(lsDomain<T, D> &passedLevelSet, std::string passedFileName)
      : levelSet(&passedLevelSet), fileName(passedFileName) {}

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  /// set file name for file to write
  void setFileName(std::string passedFileName) { fileName = passedFileName; }

  void apply() {
    // check mesh
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsReader. Not writing.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      lsMessage::getInstance()
          .addWarning("No file name specified for lsReader. Not writing.")
          .print();
      return;
    }
    if (fileName.find(".lvst") != fileName.length() - 5) {
      lsMessage::getInstance()
          .addWarning("File name does not end in '.lvst', appending it.")
          .print();
      fileName.append(".lvst");
    }

    // Open file for writing and save serialized level set in it
    std::ifstream fin(fileName);

    levelSet->deserialize(fin);

    fin.close();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsReader)

#endif // LS_READER_HPP
