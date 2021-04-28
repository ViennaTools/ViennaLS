#ifndef LS_WRITER_HPP
#define LS_WRITER_HPP

#include <fstream>

#include <lsDomain.hpp>
#include <lsPreCompileMacros.hpp>

template <class T, int D> class lsWriter {
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  std::string fileName;

public:
  lsWriter() {}

  lsWriter(lsSmartPointer<lsDomain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  lsWriter(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
           std::string passedFileName)
      : levelSet(passedLevelSet), fileName(passedFileName) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// set file name for file to write
  void setFileName(std::string passedFileName) { fileName = passedFileName; }

  void apply() {
    // check mesh
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsWriter. Not writing.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      lsMessage::getInstance()
          .addWarning("No file name specified for lsWriter. Not writing.")
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
    std::ofstream fout(fileName);

    levelSet->serialize(fout);

    fout.close();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsWriter)

#endif // LS_WRITER_HPP
