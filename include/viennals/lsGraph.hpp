#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <vcLogger.hpp>

namespace lsInternal {

using namespace viennacore;

class Graph {
  // type for indexing components
  using IndexType = std::size_t;
  // edges must be unique
  typedef std::unordered_set<IndexType> edgeListType;
  // points must be unique and addressable
  typedef std::unordered_map<IndexType, edgeListType> adjacencyListType;

  adjacencyListType adjacencyList;
  std::vector<IndexType> componentIds;

  void depthFirstComponentSearch(const std::size_t &currentVertex,
                                 int &currentComponent) {
    // set component for this vertex
    componentIds[currentVertex] = currentComponent;

    // cycle through all connected vertices and set their component
    auto vertexListIt = adjacencyList.find(currentVertex);
    if (vertexListIt == adjacencyList.end()) {
      Logger::getInstance().addError(
          "Graph: Vertex " + std::to_string(currentVertex) +
          " could not be found although it should exist!");
    }

    for (unsigned long it : vertexListIt->second) {
      // vertex has not been visited
      if (componentIds[it] == -1) {
        depthFirstComponentSearch(it, currentComponent);
      }
    }
  }

public:
  /// add new vertex
  std::size_t insertNextVertex() {
    adjacencyList.insert(std::make_pair(adjacencyList.size(), edgeListType()));
    // return new key
    return adjacencyList.size() - 1;
  }

  /// add connection of vertex1 to vertex2
  void insertNextEdge(std::size_t vertex1, std::size_t vertex2) {
    // try to insert element into vertex list
    // if it already exists, add vertex2 to vertex1 edges
    auto it =
        adjacencyList.insert(std::make_pair(vertex1, edgeListType())).first;
    it->second.insert(vertex2);

    // do the same for vertex2
    it = adjacencyList.insert(std::make_pair(vertex2, edgeListType())).first;
    it->second.insert(vertex1);
  }

  // returns a std::vector, where the value at each
  // index denotes the component, the vertex belongs to
  std::vector<IndexType> getConnectedComponents() {
    // traverse all vertices and stop at unvisited ones
    // to build connectivity.
    componentIds.resize(adjacencyList.size(), -1);
    int currentComponent = 0;
    std::size_t currentVertex = 0;

    for (auto it = adjacencyList.begin(); it != adjacencyList.end(); ++it) {
      // this vertex has not been visited, so it
      // must be part of a new component
      if (componentIds[currentVertex] == -1) {
        depthFirstComponentSearch(currentVertex, currentComponent);

        ++currentComponent;
      }
      ++currentVertex;
    }

    return componentIds;
  }

  void print() {
    std::cout << "Graph structure: " << std::endl;
    for (auto &[fst, snd] : adjacencyList) {
      auto &edges = snd;
      std::cout << "Vertex: " << fst << std::endl;
      for (const unsigned long edge : edges) {
        std::cout << edge << ", ";
      }
      std::cout << std::endl;
    }
  }
};
} // namespace lsInternal
