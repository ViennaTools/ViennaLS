#pragma once

#include <set>
#include <vector>

class lsMaterialMap {
  std::vector<int> materialMap;
  std::set<int> materials;

public:
  lsMaterialMap() = default;
  lsMaterialMap(lsMaterialMap &) = default;
  lsMaterialMap(lsMaterialMap &&) = default;

  void insertNextMaterial(const int passedMaterialId) {
    materialMap.push_back(passedMaterialId);
    materials.insert(passedMaterialId);
  }

  void setMaterialId(const std::size_t index, const int materialId) {
    if (index >= materialMap.size()) {
      materialMap.resize(index + 1);
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
};