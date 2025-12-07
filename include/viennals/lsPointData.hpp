#pragma once

#include <lsPreCompileMacros.hpp>

#include <cstdint>
#include <utility>
#include <vector>

#include <lsConcepts.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

namespace viennals {

using namespace viennacore;

/// This class holds data associated with points in space.
template <class T = double,
          lsConcepts::IsFloatingPoint<T> = lsConcepts::assignable>
class PointData {
public:
  typedef std::vector<T> ScalarDataType;
  typedef std::vector<Vec3D<T>> VectorDataType;

private:
  std::vector<ScalarDataType> scalarData;
  std::vector<std::string> scalarDataLabels;
  std::vector<VectorDataType> vectorData;
  std::vector<std::string> vectorDataLabels;

  template <class VectorType, class ReturnType = std::conditional_t<
                                  std::is_const_v<VectorType>,
                                  const typename VectorType::value_type *,
                                  typename VectorType::value_type *>>
  static ReturnType indexPointerOrNull(VectorType &v, int index) {
    if (index >= 0 && index < v.size())
      return &(v[index]);
    else
      VIENNACORE_LOG_WARNING("PointData: Tried to access out of bounds index! "
                             "Returning nullptr instead.");
    return nullptr;
  }

  template <class DataType>
  static void appendTranslateData(DataType &currentData, const DataType &source,
                                  const std::vector<unsigned> &indices) {
    currentData.reserve(currentData.size() + indices.size());
    for (unsigned int index : indices) {
      currentData.push_back(source[index]);
    }
  }

public:
  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<PointData>::New(std::forward<Args>(args)...);
  }

  /// insert new scalar data array
  void insertNextScalarData(const ScalarDataType &scalars,
                            const std::string &label = "Scalars") {
    scalarData.push_back(scalars);
    scalarDataLabels.push_back(label);
  }

  /// insert new scalar data array
  void insertNextScalarData(ScalarDataType &&scalars,
                            const std::string &label = "Scalars") {
    scalarData.push_back(std::move(scalars));
    scalarDataLabels.push_back(label);
  }

  /// insert new vector data array
  void insertNextVectorData(const VectorDataType &vectors,
                            const std::string &label = "Vectors") {
    vectorData.push_back(vectors);
    vectorDataLabels.push_back(label);
  }

  /// insert new vector data array
  void insertNextVectorData(VectorDataType &&vectors,
                            const std::string &label = "Vectors") {
    vectorData.push_back(std::move(vectors));
    vectorDataLabels.push_back(label);
  }

  /// insert or replace scalar data array
  void insertReplaceScalarData(const ScalarDataType &scalars,
                               const std::string &label = "Scalars") {
    if (int i = getScalarDataIndex(label); i != -1) {
      scalarData[i] = scalars;
    } else {
      insertNextScalarData(scalars, label);
    }
  }

  /// insert or replace scalar data array
  void insertReplaceScalarData(ScalarDataType &&scalars,
                               const std::string &label = "Scalars") {
    if (int i = getScalarDataIndex(label); i != -1) {
      scalarData[i] = std::move(scalars);
    } else {
      insertNextScalarData(std::move(scalars), label);
    }
  }

  /// insert or replace vector data array
  void insertReplaceVectorData(const VectorDataType &vectors,
                               const std::string &label = "Vectors") {
    if (int i = getVectorDataIndex(label); i != -1) {
      vectorData[i] = vectors;
    } else {
      insertNextVectorData(vectors, label);
    }
  }

  /// insert new vector data array
  void insertReplaceVectorData(VectorDataType &&vectors,
                               std::string label = "Vectors") {
    if (int i = getVectorDataIndex(label); i != -1) {
      vectorData[i] = std::move(vectors);
    } else {
      insertNextVectorData(std::move(vectors), label);
    }
  }

  /// get the number of different scalar data arrays saved
  unsigned getScalarDataSize() const { return scalarData.size(); }

  /// get the number of different vector data arrays saved
  unsigned getVectorDataSize() const { return vectorData.size(); }

  ScalarDataType *getScalarData(int index) {
    return indexPointerOrNull(scalarData, index);
  }

  const ScalarDataType *getScalarData(int index) const {
    return indexPointerOrNull(scalarData, index);
  }

  ScalarDataType *getScalarData(const std::string &searchLabel,
                                bool noWarning = false) {
    if (int i = getScalarDataIndex(searchLabel); i != -1) {
      return &(scalarData[i]);
    }
    if (!noWarning) {
      VIENNACORE_LOG_WARNING(
          "PointData attempted to access scalar data labeled '" + searchLabel +
          "', which does not exist. Returning nullptr instead.");
    }
    return nullptr;
  }

  const ScalarDataType *getScalarData(const std::string &searchLabel,
                                      bool noWarning = false) const {
    if (int i = getScalarDataIndex(searchLabel); i != -1) {
      return &(scalarData[i]);
    }
    if (!noWarning) {
      VIENNACORE_LOG_WARNING(
          "PointData attempted to access scalar data labeled '" + searchLabel +
          "', which does not exist. Returning nullptr instead.");
    }
    return nullptr;
  }

  int getScalarDataIndex(const std::string &searchLabel) const {
    for (int i = 0; i < scalarDataLabels.size(); ++i) {
      if (scalarDataLabels[i] == searchLabel) {
        return i;
      }
    }
    return -1;
  }

  std::string getScalarDataLabel(int index) const {
    return (index >= 0 && index < scalarDataLabels.size())
               ? scalarDataLabels[index]
               : "";
  }

  void setScalarDataLabel(int index, std::string newLabel) {
    scalarDataLabels[index] = std::move(newLabel);
  }

  /// Delete the scalar data at index.
  void eraseScalarData(int index) {
    scalarData.erase(scalarData.begin() + index);
    scalarDataLabels.erase(scalarDataLabels.begin() + index);
  }

  VectorDataType *getVectorData(int index) {
    return indexPointerOrNull(vectorData, index);
  }

  const VectorDataType *getVectorData(int index) const {
    return indexPointerOrNull(vectorData, index);
  }

  VectorDataType *getVectorData(const std::string &searchLabel,
                                bool noWarning = false) {
    if (int i = getVectorDataIndex(searchLabel); i != -1) {
      return &(vectorData[i]);
    }
    if (!noWarning) {
      VIENNACORE_LOG_WARNING(
          "PointData attempted to access vector data labeled '" + searchLabel +
          "', which does not exist. Returning nullptr instead.");
    }
    return nullptr;
  }

  const VectorDataType *getVectorData(const std::string &searchLabel,
                                      bool noWarning = false) const {
    if (int i = getVectorDataIndex(searchLabel); i != -1) {
      return &(vectorData[i]);
    }
    if (!noWarning) {
      VIENNACORE_LOG_WARNING(
          "PointData attempted to access vector data labeled '" + searchLabel +
          "', which does not exist. Returning nullptr instead.");
    }
    return nullptr;
  }

  int getVectorDataIndex(const std::string &searchLabel) const {
    for (int i = 0; i < vectorDataLabels.size(); ++i) {
      if (vectorDataLabels[i] == searchLabel) {
        return i;
      }
    }
    return -1;
  }

  std::string getVectorDataLabel(int index) const {
    return (index >= 0 && index < vectorDataLabels.size())
               ? vectorDataLabels[index]
               : "";
  }

  void setVectorDataLabel(int index, std::string newLabel) {
    vectorDataLabels[index] = std::move(newLabel);
  }

  /// Delete the vector data at index.
  void eraseVectorData(int index) {
    vectorData.erase(vectorData.begin() + index);
    vectorDataLabels.erase(vectorDataLabels.begin() + index);
  }

  /// Append the passed PointData to this one.
  void append(const PointData &passedData) {
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

  /// Add data in the passed source pointData into this data according to the
  /// indices passed. The index of the indices vector corresponds to the index
  /// of this data, while the values of indices correspond to the index in
  /// source.
  void translateFromData(const PointData &source,
                         const std::vector<unsigned> &indices) {
    // scalars
    for (unsigned j = 0; j < source.getScalarDataSize(); ++j) {
      insertNextScalarData(ScalarDataType(), source.getScalarDataLabel(j));
      auto currentData = --scalarData.end();
      appendTranslateData(*currentData, source.scalarData[j], indices);
    }

    // vectors
    for (unsigned j = 0; j < source.getVectorDataSize(); ++j) {
      insertNextVectorData(VectorDataType(), source.getVectorDataLabel(j));
      auto currentData = --vectorData.end();
      appendTranslateData(*currentData, source.vectorData[j], indices);
    }
  }

  /// Same as translateFromData, but the indices are given as a vector, as
  /// is the case when collecting indices during parallel algorithms.
  void translateFromMultiData(
      const PointData &source,
      const std::vector<std::vector<unsigned>> &indicesVector) {
    // scalars
    for (unsigned j = 0; j < source.getScalarDataSize(); ++j) {
      insertNextScalarData(ScalarDataType(), source.getScalarDataLabel(j));
      auto currentData = --scalarData.end();
      for (const auto &i : indicesVector) {
        appendTranslateData(*currentData, source.scalarData[j], i);
      }
    }

    // vectors
    for (unsigned j = 0; j < source.getVectorDataSize(); ++j) {
      insertNextVectorData(VectorDataType(), source.getVectorDataLabel(j));
      auto currentData = --vectorData.end();
      for (const auto &i : indicesVector) {
        appendTranslateData(*currentData, source.vectorData[j], i);
      }
    }
  }

  /// Delete all data stored in this object.
  void clear() {
    scalarData.clear();
    scalarDataLabels.clear();
    vectorData.clear();
    vectorDataLabels.clear();
  }

  /// Return whether this object is empty.
  bool empty() { return scalarData.empty() && vectorData.empty(); }

  /// Serialize PointData into a binary stream.
  std::ostream &serialize(std::ostream &stream) {
    // HEADER
    // identifier: "PointData"
    // 4 byte: number of scalar data sets
    // 4 byte: number of vector data sets
    stream << "lsPointData";
    uint32_t numberOfScalarData = scalarData.size();
    uint32_t numberOfVectorData = vectorData.size();
    stream.write(reinterpret_cast<const char *>(&numberOfScalarData),
                 sizeof(uint32_t));
    stream.write(reinterpret_cast<const char *>(&numberOfVectorData),
                 sizeof(uint32_t));

    // Scalar Data
    {
      auto labelIt = scalarDataLabels.begin();
      // iterate over all scalar data sets
      for (auto data : scalarData) {
        // write name of scalar data and size of set
        uint32_t sizeOfName = labelIt->length();
        stream.write(reinterpret_cast<char *>(&sizeOfName), sizeof(uint32_t));
        stream << *labelIt;
        uint32_t numberOfValues = data.size();
        stream.write(reinterpret_cast<const char *>(&numberOfValues),
                     sizeof(uint32_t));
        // iterate over scalars in data set
        for (auto value : data) {
          stream.write(reinterpret_cast<const char *>(&value),
                       sizeof(typename ScalarDataType::value_type));
        }
        ++labelIt;
      }
    }

    // Vector Data
    {
      auto labelIt = vectorDataLabels.begin();
      // iterate over all vector data sets
      for (auto data : vectorData) {
        // write name of vector data and size of set
        uint32_t sizeOfName = labelIt->length();
        stream.write(reinterpret_cast<char *>(&sizeOfName), sizeof(uint32_t));
        stream << *labelIt;
        uint32_t numberOfVectors = data.size();
        stream.write(reinterpret_cast<const char *>(&numberOfVectors),
                     sizeof(uint32_t));
        // iterate over vectors in data set
        for (auto vector : data) {
          // over values in vector
          for (auto value : vector) {
            stream.write(
                reinterpret_cast<const char *>(&value),
                sizeof(typename VectorDataType::value_type::value_type));
          }
        }
        ++labelIt;
      }
    }

    return stream;
  }

  /// Deserialize PointData from a binary stream.
  std::istream &deserialize(std::istream &stream) {
    char identifier[11];
    stream.read(identifier, 11);
    if (std::string(identifier).compare(0, 11, "lsPointData")) {
      Logger::getInstance()
          .addError("Reading PointData from stream failed. Header could "
                    "not be found.")
          .print();
      return stream;
    }

    // read number of data sets
    unsigned numberOfScalarData = 0;
    unsigned numberOfVectorData = 0;
    stream.read(reinterpret_cast<char *>(&numberOfScalarData),
                sizeof(uint32_t));
    stream.read(reinterpret_cast<char *>(&numberOfVectorData),
                sizeof(uint32_t));

    // read scalar data
    for (unsigned i = 0; i < numberOfScalarData; ++i) {
      uint32_t sizeOfName;
      stream.read(reinterpret_cast<char *>(&sizeOfName), sizeof(uint32_t));
      std::vector<char> dataLabel(sizeOfName);
      stream.read(dataLabel.data(), sizeOfName);
      uint32_t numberOfValues;
      stream.read(reinterpret_cast<char *>(&numberOfValues), sizeof(uint32_t));
      ScalarDataType scalarData;
      scalarData.resize(numberOfValues);
      // read all scalar values into the data
      for (auto &value : scalarData) {
        stream.read(reinterpret_cast<char *>(&value),
                    sizeof(typename ScalarDataType::value_type));
      }
      // now add this scalar data to current PointData
      insertNextScalarData(scalarData,
                           std::string(dataLabel.begin(), dataLabel.end()));
    }

    // read vector data
    for (unsigned i = 0; i < numberOfVectorData; ++i) {
      uint32_t sizeOfName;
      stream.read(reinterpret_cast<char *>(&sizeOfName), sizeof(uint32_t));
      std::vector<char> dataLabel(sizeOfName);
      stream.read(dataLabel.data(), sizeOfName);
      uint32_t numberOfValues;
      stream.read(reinterpret_cast<char *>(&numberOfValues), sizeof(uint32_t));
      VectorDataType vectorData;
      vectorData.resize(numberOfValues);
      // read all vector values into the vector data
      for (auto &vector : vectorData) {
        for (auto &value : vector) {
          stream.read(reinterpret_cast<char *>(&value),
                      sizeof(typename VectorDataType::value_type::value_type));
        }
      }
      // now add this scalar data to current PointData
      insertNextVectorData(vectorData,
                           std::string(dataLabel.begin(), dataLabel.end()));
    }

    return stream;
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION(PointData);

} // namespace viennals
