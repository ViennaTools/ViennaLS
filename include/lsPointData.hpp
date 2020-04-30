#ifndef LS_POINT_DATA_HPP
#define LS_POINT_DATA_HPP

#include <array>
#include <vector>

#include <lsMessage.hpp>

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
  /// insert new scalar data array
  void insertNextScalarData(const ScalarDataType &scalars,
                            std::string label = "Scalars") {
    scalarData.push_back(scalars);
    scalarDataLabels.push_back(label);
  }

  /// insert new vector data array
  void insertNextVectorData(const VectorDataType &vectors,
                            std::string label = "Vectors") {
    vectorData.push_back(vectors);
    vectorDataLabels.push_back(label);
  }

  /// get the number of different scalar data arrays saved
  unsigned getScalarDataSize() const { return scalarData.size(); }

  /// get the number of different vector data arrays saved
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

  bool empty() { return scalarData.empty() && vectorData.empty(); }

  /// Serialize lsPointData into a binary stream
  std::ostream &serialize(std::ostream &stream) {
    // HEADER
    // identifier: "lsPointData"
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

  std::istream &deserialize(std::istream &stream) {
    char identifier[11];
    stream.read(identifier, 11);
    if (std::string(identifier).compare(0, 11, "lsPointData")) {
      lsMessage::getInstance()
          .addWarning("Reading lsPointData from stream failed. Header could "
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
      std::vector<char> dataLabel(sizeOfName + 1);
      stream.read(dataLabel.data(), sizeOfName);
      dataLabel[sizeOfName] = '\0';
      uint32_t numberOfValues;
      stream.read(reinterpret_cast<char *>(&numberOfValues), sizeof(uint32_t));
      ScalarDataType scalarData;
      scalarData.resize(numberOfValues);
      // read all scalar values into the data
      for (auto &value : scalarData) {
        stream.read(reinterpret_cast<char *>(&value),
                    sizeof(typename ScalarDataType::value_type));
      }
      // now add this scalar data to current lsPointData
      insertNextScalarData(scalarData,
                           std::string(dataLabel.begin(), dataLabel.end()));
    }

    // read vector data
    for (unsigned i = 0; i < numberOfVectorData; ++i) {
      uint32_t sizeOfName;
      stream.read(reinterpret_cast<char *>(&sizeOfName), sizeof(uint32_t));
      std::vector<char> dataLabel(sizeOfName + 1);
      stream.read(dataLabel.data(), sizeOfName);
      dataLabel[sizeOfName] = '\0';
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
      // now add this scalar data to current lsPointData
      insertNextVectorData(vectorData,
                           std::string(dataLabel.begin(), dataLabel.end()));
    }

    return stream;
  }
};

#endif // LS_POINT_DATA_HPP
