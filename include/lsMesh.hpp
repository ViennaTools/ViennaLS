#ifndef LS_MESH_HPP
#define LS_MESH_HPP

#include <array>
#include <iostream>
#include <vector>

#include <lsPointData.hpp>

/// This class holds an explicit mesh, which is always given in 3 dimensions.
/// If it describes a 2D mesh, the third dimension is set to 0.
/// Vertices, Lines, Triangles, Tetras & Hexas are supported as geometric
/// elements.
template <class T = double> class lsMesh {
public:
  std::vector<std::array<T, 3>> nodes;
  std::vector<std::array<unsigned, 1>> vertices;
  std::vector<std::array<unsigned, 2>> lines;
  std::vector<std::array<unsigned, 3>> triangles;
  std::vector<std::array<unsigned, 4>> tetras;
  std::vector<std::array<unsigned, 8>> hexas;
  lsPointData<T> pointData;
  lsPointData<T> cellData;
  std::array<T, 3> minimumExtent;
  std::array<T, 3> maximumExtent;

private:
  // iterator typedef
  using VectorIt = typename lsPointData<T>::VectorDataType::iterator;
  // find function to avoid including the whole algorithm header
  VectorIt find(VectorIt first, VectorIt last, const std::array<T, 3> &value) {
    for (; first != last; ++first) {
      if (*first == value) {
        return first;
      }
    }
    return last;
  }

  // helper function for duplicate removal
  template <class ElementType>
  void replaceNode(ElementType &elements, std::pair<unsigned, unsigned> node) {
    for (unsigned i = 0; i < elements.size(); ++i) {
      for (unsigned j = 0; j < elements[i].size(); ++j) {
        if (elements[i][j] == node.first) {
          elements[i][j] = node.second;
        }
      }
    }
  };

public:
  const std::vector<std::array<T, 3>> &getNodes() const { return nodes; }

  std::vector<std::array<T, 3>> &getNodes() { return nodes; }

  template <int D, typename std::enable_if<D == 1, int>::type = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return vertices;
  }

  template <int D, typename std::enable_if<D == 2, int>::type = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return lines;
  }

  template <int D, typename std::enable_if<D == 3, int>::type = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return triangles;
  }

  template <int D, typename std::enable_if<D == 4, int>::type = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return tetras;
  }

  template <int D, typename std::enable_if<D == 8, int>::type = 0>
  std::vector<std::array<unsigned, D>> &getElements() {
    return hexas;
  }

  unsigned insertNextNode(const std::array<T, 3> &node) {
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
    return hexas.size();
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
    return hexas.size();
  }

  void removeDuplicateNodes() {
    std::vector<std::array<T, 3>> newNodes;
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
        duplicates.push_back(std::make_pair(i, nodeId));
      } else {
        if (adjusted)
          duplicates.push_back(std::make_pair(i, newNodes.size()));
        newNodes.push_back(nodes[i]);
      }
    }
    nodes = newNodes;

    // now replace in vertices
    // TODO also need to shift down all other nodes
    for (unsigned i = 0; i < duplicates.size(); ++i) {
      replaceNode(vertices, duplicates[i]);
      replaceNode(lines, duplicates[i]);
      replaceNode(triangles, duplicates[i]);
      replaceNode(tetras, duplicates[i]);
      replaceNode(hexas, duplicates[i]);
    }
  }

  void append(const lsMesh<T> &passedMesh) {
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
    std::cout << "lsMesh:" << std::endl;
    std::cout << "Number of Nodes: " << nodes.size() << std::endl;
    if (vertices.size() > 0)
      std::cout << "Number of Vertices: " << vertices.size() << std::endl;
    if (lines.size() > 0)
      std::cout << "Number of Lines: " << lines.size() << std::endl;
    if (triangles.size() > 0)
      std::cout << "Number of Triangles: " << triangles.size() << std::endl;
    if (tetras.size() > 0)
      std::cout << "Number of Tetrahedrons: " << tetras.size() << std::endl;
    if (hexas.size() > 0)
      std::cout << "Number of Hexas: " << hexas.size() << std::endl;
    // pointData
    if (pointData.getScalarDataSize() > 0) {
      std::cout << "Scalar data:" << std::endl;
      for (unsigned i = 0; i < pointData.getScalarDataSize(); ++i) {
        std::cout << "  \"" << pointData.getScalarDataLabel(i)
                  << "\" of size " << pointData.getScalarData(i)->size()
                  << std::endl;
      }
    }
    if (pointData.getVectorDataSize() > 0) {
      std::cout << "Vector data:" << std::endl;
      for (unsigned i = 0; i < pointData.getVectorDataSize(); ++i) {
        std::cout << "  \"" << pointData.getVectorDataLabel(i)
                  << "\" of size " << pointData.getVectorData(i)->size()
                  << std::endl;
      }
    }


    // cellData
    if (cellData.getScalarDataSize() > 0) {
      std::cout << "Scalar data:" << std::endl;
      for (unsigned i = 0; i < cellData.getScalarDataSize(); ++i) {
        std::cout << "  \"" << cellData.getScalarDataLabel(i)
                  << "\" of size " << cellData.getScalarData(i)->size()
                  << std::endl;
      }
    }
    if (cellData.getVectorDataSize() > 0) {
      std::cout << "Vector data:" << std::endl;
      for (unsigned i = 0; i < cellData.getVectorDataSize(); ++i) {
        std::cout << "  \"" << cellData.getVectorDataLabel(i)
                  << "\" of size " << cellData.getVectorData(i)->size()
                  << std::endl;
      }
    }
  }
};

#endif // LS_MESH_HPP

// add all template specialisations for this class
PRECOMPILE_PRECISION(lsMesh);