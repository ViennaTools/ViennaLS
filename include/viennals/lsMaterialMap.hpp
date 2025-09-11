#pragma once

#include <set>
#include <vector>

namespace viennals {

/// \brief A class for mapping layer indices to material IDs
///
/// This class maintains a mapping between layer indices and material IDs,
/// keeping track of both the sequential material assignment and the unique
/// materials used in the mapping.
class MaterialMap {
  std::vector<int> materialMap;
  std::set<int> materials;

public:
  MaterialMap() = default;
  MaterialMap(const MaterialMap &) = default;
  MaterialMap &operator=(const MaterialMap &) = default;
  MaterialMap &operator=(MaterialMap &&) = default;

  void insertNextMaterial(const int passedMaterialId) {
    materialMap.push_back(passedMaterialId);
    materials.insert(passedMaterialId);
  }

  void setMaterialId(const std::size_t index, const int materialId) {
    if (index >= materialMap.size()) {
      materialMap.resize(index + 1, -1); // Initialize new elements with -1
    }
    materialMap[index] = materialId;
    materials.insert(materialId);
  }

  std::size_t getNumberOfLayers() const { return materialMap.size(); }

  std::size_t getNumberOfMaterials() const { return materials.size(); }

  int getMaterialId(const std::size_t index) const {
    if (index >= materialMap.size())
      return -1;
    return materialMap[index];
  }

  // Additional utility methods
  bool isValidIndex(const std::size_t index) const {
    return index < materialMap.size();
  }

  void clear() {
    materialMap.clear();
    materials.clear();
  }

  void reserve(const std::size_t size) { materialMap.reserve(size); }

  bool hasMaterial(const int materialId) const {
    return materials.count(materialId) > 0;
  }

  const std::set<int> &getMaterials() const { return materials; }

  const std::vector<int> &getMaterialMap() const { return materialMap; }
};

} // namespace viennals
