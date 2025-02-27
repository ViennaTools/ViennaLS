#pragma once

#include <lsPreCompileMacros.hpp>

#include <array>
#include <iostream>
#include <vector>

#include <lsPointData.hpp>

#include <vcVectorUtil.hpp>

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
  Vec3D<T> minimumExtent;
  Vec3D<T> maximumExtent;

private:
  // iterator typedef
  using VectorIt = typename PointData<T>::VectorDataType::iterator;
  // find function to avoid including the whole algorithm header
  static VectorIt find(VectorIt first, VectorIt last, const Vec3D<T> &value) {
    for (; first != last; ++first) {
      if (*first == value) {
        return first;
      }
    }
    return last;
  }

  // helper function for duplicate removal
  template <class ElementType>
  static void replaceNode(ElementType &elements,
                          std::pair<unsigned, unsigned> node) {
    for (unsigned i = 0; i < elements.size(); ++i) {
      for (unsigned j = 0; j < elements[i].size(); ++j) {
        if (elements[i][j] == node.first) {
          elements[i][j] = node.second;
        }
      }
    }
  };

public:
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

  void removeDuplicateNodes() {
    std::vector<Vec3D<T>> newNodes;
    // can just push first point since it cannot be duplicate
    newNodes.push_back(nodes[0]);
    // now check for duplicates
    // pair of oldId <-> newId
    std::vector<std::pair<unsigned, unsigned>> duplicates;
    bool adjusted = false;
    for (unsigned i = 1; i < nodes.size(); ++i) {
      auto it = find(newNodes.begin(), newNodes.end(), nodes[i]);
      if (it != newNodes.end()) {
        adjusted = true;
        // if duplicate point, save it to be replaced
        unsigned nodeId = std::distance(newNodes.begin(), it);
        duplicates.emplace_back(i, nodeId);
      } else {
        if (adjusted)
          duplicates.push_back(std::make_pair(i, newNodes.size()));
        newNodes.push_back(nodes[i]);
      }
    }
    nodes = newNodes;

    // now replace in vertices
    // TODO also need to shift down all other nodes
    for (auto &duplicate : duplicates) {
      replaceNode(vertices, duplicate);
      replaceNode(lines, duplicate);
      replaceNode(triangles, duplicate);
      replaceNode(tetras, duplicate);
      replaceNode(hexas, duplicate);
    }
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