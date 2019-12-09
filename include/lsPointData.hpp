#ifndef LS_POINT_DATA_HPP
#define LS_POINT_DATA_HPP

#include <array>
#include <vector>

/// This class holds data associated with points in space.
class lsPointData {
public:
  typedef std::vector<double> ScalarDataType;
  typedef std::vector<std::array<double, 3>> VectorDataType;

private:
  std::vector<ScalarDataType> scalarData;
  std::vector<std::string> scalarDataLabels;
  std::vector<VectorDataType> vectorData;
  std::vector<std::string> vectorDataLabels;

public:
  void insertNextScalarData(const ScalarDataType &scalars,
                            std::string label = "Scalars") {
    scalarData.push_back(scalars);
    scalarDataLabels.push_back(label);
  }

  void insertNextVectorData(const VectorDataType &vectors,
                            std::string label = "Vectors") {
    vectorData.push_back(vectors);
    vectorDataLabels.push_back(label);
  }

  unsigned getScalarDataSize() const { return scalarData.size(); }

  unsigned getVectorDataSize() const { return vectorData.size(); }

  ScalarDataType *getScalarData(int index) { return &(scalarData[index]); }

  const ScalarDataType *getScalarData(int index) const {
    return &(scalarData[index]);
  }

  ScalarDataType *getScalarData(std::string searchLabel) {
    for (unsigned i = 0; i < scalarDataLabels.size(); ++i) {
      if (scalarDataLabels[i] == searchLabel) {
        return &(scalarData[i]);
      }
    }
    return nullptr;
  }

  const ScalarDataType *getScalarData(std::string searchLabel) const {
    for (unsigned i = 0; i < scalarDataLabels.size(); ++i) {
      if (scalarDataLabels[i] == searchLabel) {
        return &(scalarData[i]);
      }
    }
    return nullptr;
  }

  std::string getScalarDataLabel(int index) const {
    return scalarDataLabels[index];
  }

  VectorDataType *getVectorData(int index) { return &(vectorData[index]); }

  const VectorDataType *getVectorData(int index) const {
    return &(vectorData[index]);
  }

  VectorDataType *getVectorData(std::string searchLabel) {
    for (unsigned i = 0; i < vectorDataLabels.size(); ++i) {
      if (vectorDataLabels[i] == searchLabel) {
        return &(vectorData[i]);
      }
    }
    return nullptr;
  }

  const VectorDataType *getVectorData(std::string searchLabel) const {
    for (unsigned i = 0; i < vectorDataLabels.size(); ++i) {
      if (vectorDataLabels[i] == searchLabel) {
        return &(vectorData[i]);
      }
    }
    return nullptr;
  }

  std::string getVectorDataLabel(int index) const {
    return vectorDataLabels[index];
  }

  void append(const lsPointData &passedData) {
    scalarData.insert(scalarData.end(), passedData.scalarData.begin(),
                      passedData.scalarData.end());
    scalarDataLabels.insert(scalarDataLabels.end(),
                            passedData.scalarDataLabels.begin(),
                            passedData.scalarDataLabels.end());
    vectorData.insert(vectorData.end(), passedData.vectorData.begin(),
                      passedData.vectorData.end());
    vectorDataLabels.insert(vectorDataLabels.end(),
                            passedData.vectorDataLabels.begin(),
                            passedData.vectorDataLabels.end());
  }

  void clear() {
    scalarData.clear();
    scalarDataLabels.clear();
    vectorData.clear();
    vectorDataLabels.clear();
  }
};

#endif // LS_POINT_DATA_HPP
