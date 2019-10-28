#ifndef LS_GRAPH_HPP
#define LS_GRAPH_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <lsMessage.hpp>

namespace lsInternal {

class lsGraph {
  // edges must be unique
  typedef typename std::unordered_set<std::size_t> edgeListType;
  // points must be unique and adressable
  typedef typename std::unordered_map<std::size_t, edgeListType>
      adjacencyListType;

  adjacencyListType adjacencyList;
  std::vector<int> componentIds;

  void depthFirstComponentSearch(const std::size_t &currentVertex,
                                 int &currentComponent) {
    // set component for this vertex
    componentIds[currentVertex] = currentComponent;

    // cylce through all connected vertices and set their component
    auto vertexListIt = adjacencyList.find(currentVertex);
    if (vertexListIt == adjacencyList.end()) {
      lsMessage::getInstance().addError(
          "lsGraph: Vertex " + std::to_string(currentVertex) +
          " could not be found although it should exist!");
    }

    for (auto it = vertexListIt->second.begin();
         it != vertexListIt->second.end(); ++it) {
      // vertex has not been visited
      if (componentIds[*it] == -1) {
        depthFirstComponentSearch(*it, currentComponent);
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

  // returns an std::vector, where the value at each
  // index denotes the component, the vertex belongs to
  std::vector<int> getConnectedComponents() {
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
    for (auto it = adjacencyList.begin(); it != adjacencyList.end(); ++it) {
      auto &edges = (*it).second;
      std::cout << "Vertex: " << (*it).first << std::endl;
      for (auto edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt) {
        std::cout << *edgeIt << ", ";
      }
      std::cout << std::endl;
    }
  }
};
} // namespace lsInternal

#endif // LS_GRAPH_HPP
