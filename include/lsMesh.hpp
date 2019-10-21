#ifndef LS_MESH_HPP
#define LS_MESH_HPP

#include <iostream>
#include <vector>

#include <hrleVectorType.hpp>

/// This class holds an explicit mesh, which is always given in 3 dimensions.
/// If it describes a 2D mesh, the third dimension is set to 0.
/// Vertices, Lines, Triangles, Tetras & Hexas are supported as geometric
/// elements.
class lsMesh {
public:
  std::vector<hrleVectorType<double, 3>> nodes;
  std::vector<hrleVectorType<unsigned, 1>> vertices;
  std::vector<hrleVectorType<unsigned, 2>> lines;
  std::vector<hrleVectorType<unsigned, 3>> triangles;
  std::vector<hrleVectorType<unsigned, 4>> tetras;
  std::vector<hrleVectorType<unsigned, 8>> hexas;
  std::vector<unsigned> materials;
  std::vector<std::vector<double>> scalarData;
  std::vector<std::string> scalarDataLabels;
  std::vector<std::vector<hrleVectorType<double, 3>>> vectorData;
  std::vector<std::string> vectorDataLabels;

private:
  // helper function for duplicate removal
  template <class ElementType>
  void replaceNode(ElementType &elements, std::pair<unsigned, unsigned> node) {
    for (unsigned i = 0; i < elements.size(); ++i) {
      for (unsigned j = 0; j < ElementType::value_type::dimension; ++j) {
        if (elements[i][j] == node.first) {
          elements[i][j] = node.second;
        }
      }
    }
  };

public:
  const std::vector<hrleVectorType<double, 3>> &getNodes() const {
    return nodes;
  }

  std::vector<hrleVectorType<double, 3>> &getNodes() { return nodes; }

  std::vector<double> &getScalarData(int i) { return scalarData[i]; }

  std::vector<hrleVectorType<double, 3>> &getVectorData(int i) {
    return vectorData[i];
  }

  std::vector<unsigned> &getMaterials() { return materials; }

  template <int D, typename std::enable_if<D == 1, int>::type = 0>
  std::vector<hrleVectorType<unsigned, D>> &getElements() {
    return vertices;
  }

  template <int D, typename std::enable_if<D == 2, int>::type = 0>
  std::vector<hrleVectorType<unsigned, D>> &getElements() {
    return lines;
  }

  template <int D, typename std::enable_if<D == 3, int>::type = 0>
  std::vector<hrleVectorType<unsigned, D>> &getElements() {
    return triangles;
  }

  template <int D, typename std::enable_if<D == 4, int>::type = 0>
  std::vector<hrleVectorType<unsigned, D>> &getElements() {
    return tetras;
  }

  template <int D, typename std::enable_if<D == 8, int>::type = 0>
  std::vector<hrleVectorType<unsigned, D>> &getElements() {
    return hexas;
  }

  unsigned insertNextNode(hrleVectorType<double, 3> &node) {
    nodes.push_back(node);
    return nodes.size() - 1;
  }

  unsigned insertNextVertex(hrleVectorType<unsigned, 1> &vertex) {
    vertices.push_back(vertex);
    return vertices.size() - 1;
  }

  unsigned insertNextLine(hrleVectorType<unsigned, 2> &line) {
    lines.push_back(line);
    return lines.size() - 1;
  }

  unsigned insertNextTriangle(hrleVectorType<unsigned, 3> &triangle) {
    triangles.push_back(triangle);
    return triangles.size() - 1;
  }

  unsigned insertNextTetra(hrleVectorType<unsigned, 4> &tetra) {
    tetras.push_back(tetra);
    return tetras.size() - 1;
  }

  unsigned insertNextHexa(hrleVectorType<unsigned, 8> &hexa) {
    hexas.push_back(hexa);
    return hexas.size();
  }

  unsigned insertNextElement(hrleVectorType<unsigned, 1> &vertex) {
    vertices.push_back(vertex);
    return vertices.size() - 1;
  }

  unsigned insertNextElement(hrleVectorType<unsigned, 2> &line) {
    lines.push_back(line);
    return lines.size() - 1;
  }

  unsigned insertNextElement(hrleVectorType<unsigned, 3> &triangle) {
    triangles.push_back(triangle);
    return triangles.size() - 1;
  }

  unsigned insertNextElement(hrleVectorType<unsigned, 4> &tetra) {
    tetras.push_back(tetra);
    return tetras.size() - 1;
  }

  unsigned insertNextElement(hrleVectorType<unsigned, 8> &hexa) {
    hexas.push_back(hexa);
    return hexas.size();
  }

  void insertNextScalarData(std::vector<double> &scalars,
                            std::string label = "Scalars") {
    scalarData.push_back(scalars);
    scalarDataLabels.push_back(label);
  }

  void insertNextVectorData(std::vector<hrleVectorType<double, 3>> &vectors,
                            std::string label = "Vectors") {
    vectorData.push_back(vectors);
    vectorDataLabels.push_back(label);
  }

  void removeDuplicateNodes() {
    std::vector<hrleVectorType<double, 3>> newNodes;
    // can just push first point since it cannot be duplicate
    newNodes.push_back(nodes[0]);
    // now check for duplicates
    // pair of oldId <-> newId
    std::vector<std::pair<unsigned, unsigned>> duplicates;
    bool adjusted = false;
    for (unsigned i = 1; i < nodes.size(); ++i) {
      auto it = std::find(newNodes.begin(), newNodes.end(), nodes[i]);
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

  void clear() {
    nodes.clear();
    vertices.clear();
    lines.clear();
    triangles.clear();
    tetras.clear();
    hexas.clear();
    materials.clear();
    scalarData.clear();
    scalarDataLabels.clear();
    vectorData.clear();
    vectorDataLabels.clear();
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
  }
};

#endif // LS_MESH_HPP
