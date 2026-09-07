#pragma once

#include <lsPreCompileMacros.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <vcPointData.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// This class holds an explicit mesh, which is always given in 3 dimensions.
/// If it describes a 2D mesh, the third dimension is set to 0.
/// Vertices, Lines, Triangles, Tetras & Hexas are supported as geometric
/// elements.
template <class T = double> class Mesh {
public:
  std::vector<Vec3D<T>> nodes;
  std::vector<std::array<unsigned, 1>> vertices;
  std::vector<std::array<unsigned, 2>> lines;
  std::vector<std::array<unsigned, 3>> triangles;
  std::vector<std::array<unsigned, 4>> tetras;
  std::vector<std::array<unsigned, 8>> hexas;
  PointData<T> pointData;
  PointData<T> cellData;
  Vec3D<T> minimumExtent{};
  Vec3D<T> maximumExtent{};

  constexpr static const char *materialIdsLabel = "MaterialIds";
  constexpr static const char *normalsLabel = "Normals";

  // Convenience function to create a new mesh smart pointer.
  static auto New() { return SmartPointer<Mesh>::New(); }

  const std::vector<Vec3D<T>> &getNodes() const { return nodes; }

  std::vector<Vec3D<T>> &getNodes() { return nodes; }

  template <int D, std::enable_if_t<D == 1, int> = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return vertices;
  }

  template <int D, std::enable_if_t<D == 2, int> = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return lines;
  }

  template <int D, std::enable_if_t<D == 3, int> = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return triangles;
  }

  template <int D, std::enable_if_t<D == 4, int> = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return tetras;
  }

  template <int D, std::enable_if_t<D == 8, int> = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return hexas;
  }

  PointData<T> &getPointData() { return pointData; }

  const PointData<T> &getPointData() const { return pointData; }

  PointData<T> &getCellData() { return cellData; }

  const PointData<T> &getCellData() const { return cellData; }

  // helper function to get normals
  std::vector<Vec3D<T>> *getNormals(const char *label = normalsLabel) {
    return cellData.getVectorData(label);
  }

  // helper function to get material ids
  std::vector<T> *getMaterialIds(const char *label = materialIdsLabel) {
    return cellData.getScalarData(label);
  }

  unsigned insertNextNode(const Vec3D<T> &node) {
    nodes.push_back(node);
    return nodes.size() - 1;
  }

  unsigned insertNextVertex(const std::array<unsigned, 1> &vertex) {
    vertices.push_back(vertex);
    return vertices.size() - 1;
  }

  unsigned insertNextLine(const std::array<unsigned, 2> &line) {
    lines.push_back(line);
    return lines.size() - 1;
  }

  unsigned insertNextTriangle(const std::array<unsigned, 3> &triangle) {
    triangles.push_back(triangle);
    return triangles.size() - 1;
  }

  unsigned insertNextTetra(const std::array<unsigned, 4> &tetra) {
    tetras.push_back(tetra);
    return tetras.size() - 1;
  }

  unsigned insertNextHexa(const std::array<unsigned, 8> &hexa) {
    hexas.push_back(hexa);
    return hexas.size() - 1;
  }

  unsigned insertNextElement(const std::array<unsigned, 1> &vertex) {
    vertices.push_back(vertex);
    return vertices.size() - 1;
  }

  unsigned insertNextElement(const std::array<unsigned, 2> &line) {
    lines.push_back(line);
    return lines.size() - 1;
  }

  unsigned insertNextElement(const std::array<unsigned, 3> &triangle) {
    triangles.push_back(triangle);
    return triangles.size() - 1;
  }

  unsigned insertNextElement(const std::array<unsigned, 4> &tetra) {
    tetras.push_back(tetra);
    return tetras.size() - 1;
  }

  unsigned insertNextElement(const std::array<unsigned, 8> &hexa) {
    hexas.push_back(hexa);
    return hexas.size() - 1;
  }

  /// Remove exactly equal nodes, preserving first-occurrence order and the
  /// first node's scalar/vector point data. Remap all element node IDs;
  /// elements and cell data are retained, including degenerate elements.
  /// Nodes containing NaNs remain distinct. Expected linear time in the number
  /// of nodes, connectivity entries, and point-data values, with linear scratch
  /// storage. When nodes are merged, point-data arrays must contain one value
  /// per input node; otherwise throws std::invalid_argument without changing
  /// the mesh.
  void removeDuplicateNodes() {
    if (nodes.size() < 2)
      return;

    struct NodeHash {
      std::size_t operator()(const Vec3D<T> &node) const {
        std::size_t seed = 0;
        for (const T coordinate : node) {
          seed ^= std::hash<T>{}(coordinate) + std::size_t(0x9e3779b9) +
                  (seed << 6) + (seed >> 2);
        }
        return seed;
      }
    };

    // Keep full coordinates as keys so hash collisions cannot merge nodes.
    // std::hash<T> also gives equal hashes for +0 and -0, which compare equal.
    std::unordered_map<Vec3D<T>, unsigned, NodeHash> uniqueNodes;
    uniqueNodes.reserve(nodes.size());
    std::vector<Vec3D<T>> newNodes;
    newNodes.reserve(nodes.size());
    std::vector<unsigned> oldToNew(nodes.size());
    std::vector<unsigned> retainedIndices;
    retainedIndices.reserve(nodes.size());

    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const auto &node = nodes[i];
      const auto newId = static_cast<unsigned>(newNodes.size());
      // NaNs are not equal to themselves and cannot serve as hash-table keys.
      if (!std::isnan(node[0]) && !std::isnan(node[1]) &&
          !std::isnan(node[2])) {
        const auto result = uniqueNodes.try_emplace(node, newId);
        oldToNew[i] = result.first->second;
        if (!result.second)
          continue;
      } else {
        oldToNew[i] = newId;
      }
      newNodes.push_back(node);
      retainedIndices.push_back(static_cast<unsigned>(i));
    }

    if (newNodes.size() == nodes.size())
      return;

    // Validate before translating data or changing any mesh connectivity.
    const auto validateData = [this](const auto &arrays) {
      for (const auto &data : arrays) {
        if (data.size() != nodes.size()) {
          throw std::invalid_argument(
              "Mesh::removeDuplicateNodes: point-data size must match the "
              "number of nodes.");
        }
      }
    };
    validateData(pointData.getScalarData());
    validateData(pointData.getVectorData());
    PointData<T> newPointData;
    newPointData.translateFromData(pointData, retainedIndices);

    const auto remapElements = [&oldToNew](auto &elements) {
      for (auto &element : elements) {
        for (auto &nodeId : element)
          nodeId = oldToNew[nodeId];
      }
    };
    remapElements(vertices);
    remapElements(lines);
    remapElements(triangles);
    remapElements(tetras);
    remapElements(hexas);

    nodes = std::move(newNodes);
    pointData = std::move(newPointData);
  }

  void append(const Mesh<T> &passedMesh) {
    const unsigned numberOfOldNodes = nodes.size();

    // append new nodes
    nodes.insert(nodes.end(), passedMesh.nodes.begin(), passedMesh.nodes.end());

    // go through all elements and increase node IDs to match new IDS
    const unsigned numberOfVertices = vertices.size();
    vertices.insert(vertices.end(), passedMesh.vertices.begin(),
                    passedMesh.vertices.end());
    for (unsigned i = numberOfVertices;
         i < passedMesh.vertices.size() + numberOfVertices; ++i) {
      vertices[i][0] += numberOfOldNodes;
    }

    const unsigned numberOfLines = lines.size();
    lines.insert(lines.end(), passedMesh.lines.begin(), passedMesh.lines.end());
    for (unsigned i = numberOfLines;
         i < passedMesh.lines.size() + numberOfLines; ++i) {
      for (unsigned d = 0; d < 2; ++d) {
        lines[i][d] += numberOfOldNodes;
      }
    }

    const unsigned numberOfTriangles = triangles.size();
    triangles.insert(triangles.end(), passedMesh.triangles.begin(),
                     passedMesh.triangles.end());
    for (unsigned i = numberOfTriangles;
         i < passedMesh.triangles.size() + numberOfTriangles; ++i) {
      for (unsigned d = 0; d < 3; ++d) {
        triangles[i][d] += numberOfOldNodes;
      }
    }

    const unsigned numberOfTetras = tetras.size();
    tetras.insert(tetras.end(), passedMesh.tetras.begin(),
                  passedMesh.tetras.end());
    for (unsigned i = numberOfTetras;
         i < passedMesh.tetras.size() + numberOfTetras; ++i) {
      for (unsigned d = 0; d < 4; ++d) {
        tetras[i][d] += numberOfOldNodes;
      }
    }

    const unsigned numberOfHexas = hexas.size();
    hexas.insert(hexas.end(), passedMesh.hexas.begin(), passedMesh.hexas.end());
    for (unsigned i = numberOfHexas;
         i < passedMesh.hexas.size() + numberOfHexas; ++i) {
      for (unsigned d = 0; d < 8; ++d) {
        hexas[i][d] += numberOfOldNodes;
      }
    }

    // Append data
    // TODO need to adjust lsVTKWriter to deal with different data correctly
    // currently this only works for vertex only meshes
    pointData.append(passedMesh.pointData);
    cellData.append(passedMesh.cellData);

    // if(lsPointData<T>::scalarData.size() < nodes.size())
    for (unsigned i = 0; i < pointData.getScalarDataSize(); ++i) {
      pointData.getScalarData(i)->resize(vertices.size());
    }
    for (unsigned i = 0; i < pointData.getVectorDataSize(); ++i) {
      pointData.getVectorData(i)->resize(vertices.size());
    }

    for (unsigned i = 0; i < cellData.getScalarDataSize(); ++i) {
      cellData.getScalarData(i)->resize(vertices.size());
    }
    for (unsigned i = 0; i < cellData.getVectorDataSize(); ++i) {
      cellData.getVectorData(i)->resize(vertices.size());
    }
  }

  void clear() {
    nodes.clear();
    vertices.clear();
    lines.clear();
    triangles.clear();
    tetras.clear();
    hexas.clear();
    pointData.clear();
    cellData.clear();
    minimumExtent = Vec3D<T>{};
    maximumExtent = Vec3D<T>{};
  }

  void print() {
    std::cout << "Mesh:" << std::endl;
    std::cout << "Number of Nodes: " << nodes.size() << std::endl;
    if (!vertices.empty())
      std::cout << "Number of Vertices: " << vertices.size() << std::endl;
    if (!lines.empty())
      std::cout << "Number of Lines: " << lines.size() << std::endl;
    if (!triangles.empty())
      std::cout << "Number of Triangles: " << triangles.size() << std::endl;
    if (!tetras.empty())
      std::cout << "Number of Tetrahedrons: " << tetras.size() << std::endl;
    if (!hexas.empty())
      std::cout << "Number of Hexas: " << hexas.size() << std::endl;
    // pointData
    if (pointData.getScalarDataSize() > 0) {
      std::cout << "Scalar data:" << std::endl;
      for (unsigned i = 0; i < pointData.getScalarDataSize(); ++i) {
        std::cout << "  \"" << pointData.getScalarDataLabel(i) << "\" of size "
                  << pointData.getScalarData(i)->size() << std::endl;
      }
    }
    if (pointData.getVectorDataSize() > 0) {
      std::cout << "Vector data:" << std::endl;
      for (unsigned i = 0; i < pointData.getVectorDataSize(); ++i) {
        std::cout << "  \"" << pointData.getVectorDataLabel(i) << "\" of size "
                  << pointData.getVectorData(i)->size() << std::endl;
      }
    }

    // cellData
    if (cellData.getScalarDataSize() > 0) {
      std::cout << "Scalar data:" << std::endl;
      for (unsigned i = 0; i < cellData.getScalarDataSize(); ++i) {
        std::cout << "  \"" << cellData.getScalarDataLabel(i) << "\" of size "
                  << cellData.getScalarData(i)->size() << std::endl;
      }
    }
    if (cellData.getVectorDataSize() > 0) {
      std::cout << "Vector data:" << std::endl;
      for (unsigned i = 0; i < cellData.getVectorDataSize(); ++i) {
        std::cout << "  \"" << cellData.getVectorDataLabel(i) << "\" of size "
                  << cellData.getVectorData(i)->size() << std::endl;
      }
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION(Mesh);

} // namespace viennals
