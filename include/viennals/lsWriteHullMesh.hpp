#pragma once

#ifdef VIENNALS_USE_VTK // this class needs vtk support

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include <hrleDenseIterator.hpp>

#include <lsDomain.hpp>
#include <lsMaterialMap.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <unordered_map>
#include <utility>

namespace viennals {

using namespace viennacore;

/// Extracts a hull (surface outline) mesh from a stack of level sets using
/// ToMultiSurfaceMesh with closed boundary caps. Material IDs are assigned per
/// cell. setSharpCorners(true/false) controls whether sharp corner generation
/// is enabled during surface extraction. The result is written to
/// fileName + "_hull.vtp".
template <class T, int D> class WriteHullMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;
  using LevelSetsType = std::vector<SmartPointer<Domain<T, D>>>;
  LevelSetsType levelSets;
  SmartPointer<MaterialMap> materialMap = nullptr;
  std::string fileName;
  bool writeToFile = true;
  bool generateSharpCorners = false;
  std::unordered_map<std::string, std::vector<double>> metaData;

  vtkSmartPointer<vtkPolyData> hullVTK;

  void addMetaDataToVTK(vtkDataSet *data) const {
    if (metaData.empty())
      return;

    vtkSmartPointer<vtkFieldData> fieldData = data->GetFieldData();
    for (const auto &meta : metaData) {
      if (meta.second.empty())
        continue;

      vtkSmartPointer<vtkFloatArray> metaDataArray =
          vtkSmartPointer<vtkFloatArray>::New();
      metaDataArray->SetName(meta.first.c_str());
      metaDataArray->SetNumberOfValues(meta.second.size());
      for (size_t i = 0; i < meta.second.size(); ++i)
        metaDataArray->SetValue(i, meta.second[i]);
      fieldData->AddArray(metaDataArray);
    }
  }

  /// Generate hull mesh using ToMultiSurfaceMesh and boundary caps.
  /// generateSharpCorners is forwarded to the converter. Works for both
  /// 2D (line segments) and 3D (triangles).
  void generateSharpHull() {
    auto &grid = levelSets[0]->getGrid();
    double gridDelta = grid.getGridDelta();

    int totalMinimum[D];
    int totalMaximum[D];
    for (int i = 0; i < D; ++i) {
      totalMinimum[i] = std::numeric_limits<int>::max();
      totalMaximum[i] = std::numeric_limits<int>::lowest();
    }

    for (auto &ls : levelSets) {
      if (ls->getNumberOfPoints() == 0)
        continue;
      auto &dom = ls->getDomain();
      auto &grd = ls->getGrid();
      for (int i = 0; i < D; ++i) {
        int minP = (grd.isNegBoundaryInfinite(i)) ? dom.getMinRunBreak(i) - 1
                                                  : grd.getMinBounds(i);
        int maxP = (grd.isPosBoundaryInfinite(i)) ? dom.getMaxRunBreak(i) + 1
                                                  : grd.getMaxBounds(i);
        if (minP < totalMinimum[i])
          totalMinimum[i] = minP;
        if (maxP > totalMaximum[i])
          totalMaximum[i] = maxP;
      }
    }

    // Generate multi-surface mesh with sharp corners
    auto multiMesh = SmartPointer<Mesh<T>>::New();
    ToMultiSurfaceMesh<T, D> converter;
    for (auto &ls : levelSets)
      converter.insertNextLevelSet(ls);
    if (materialMap)
      converter.setMaterialMap(materialMap);
    converter.setMesh(multiMesh);
    converter.setSharpCorners(generateSharpCorners);
    converter.apply();

    // Add closed boundary caps
    if constexpr (D == 2) {
      auto capBoundary = [&](int boundaryAxis, double fixedValue,
                             double rangeMin, double rangeMax) {
        int varyingAxis = 1 - boundaryAxis;
        std::vector<std::pair<double, unsigned>> boundaryNodes;
        double tolerance = 1e-6 * gridDelta;

        for (unsigned i = 0; i < multiMesh->nodes.size(); ++i) {
          if (std::abs(multiMesh->nodes[i][boundaryAxis] - fixedValue) <
              tolerance) {
            if (multiMesh->nodes[i][varyingAxis] >= rangeMin - tolerance &&
                multiMesh->nodes[i][varyingAxis] <= rangeMax + tolerance) {
              boundaryNodes.push_back({multiMesh->nodes[i][varyingAxis], i});
            }
          }
        }

        auto addCorner = [&](double val) {
          for (const auto &bn : boundaryNodes)
            if (std::abs(bn.first - val) < tolerance)
              return;

          if (levelSets.back()->getNumberOfPoints() == 0)
            return;

          viennahrle::Index<D> idx;
          idx[boundaryAxis] = std::round(fixedValue / gridDelta);
          idx[varyingAxis] = std::round(val / gridDelta);

          auto &grid = levelSets.back()->getGrid();
          if (grid.isOutsideOfDomain(idx))
            idx = grid.globalIndices2LocalIndices(idx);

          // Clamp to grid bounds to ensure iterator safety
          auto minGrid = grid.getMinGridPoint();
          auto maxGrid = grid.getMaxGridPoint();
          for (int i = 0; i < D; ++i) {
            if (idx[i] < minGrid[i])
              idx[i] = minGrid[i];
            if (idx[i] > maxGrid[i])
              idx[i] = maxGrid[i];
          }

          viennahrle::ConstDenseIterator<hrleDomainType> it(
              levelSets.back()->getDomain());
          it.goToIndices(idx);
          if (it.getValue() > 0.)
            return;

          Vec3D<T> pos;
          pos[boundaryAxis] = fixedValue;
          pos[varyingAxis] = val;
          pos[2] = 0;
          unsigned id = multiMesh->insertNextNode(pos);
          boundaryNodes.push_back({val, id});
        };
        addCorner(rangeMin);
        addCorner(rangeMax);

        std::sort(boundaryNodes.begin(), boundaryNodes.end());

        if (boundaryNodes.size() < 2)
          return;

        auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
        for (size_t i = 0; i < boundaryNodes.size() - 1; ++i) {
          double mid =
              (boundaryNodes[i].first + boundaryNodes[i + 1].first) * 0.5;

          int matId = -1;
          int boundaryCandidate = -1;
          viennahrle::Index<D> queryIdx;
          queryIdx[boundaryAxis] = std::round(fixedValue / gridDelta);
          queryIdx[varyingAxis] = std::round(mid / gridDelta);

          for (unsigned l = 0; l < levelSets.size(); ++l) {
            viennahrle::ConstDenseIterator<hrleDomainType> it(
                levelSets[l]->getDomain());
            it.goToIndices(queryIdx);
            double val = it.getValue();
            if (val < -1e-6 * gridDelta) {
              matId = (materialMap) ? materialMap->getMaterialId(l) : (int)l;
              break;
            }
            if (val <= 1e-6 * gridDelta && boundaryCandidate == -1) {
              boundaryCandidate =
                  (materialMap) ? materialMap->getMaterialId(l) : (int)l;
            }
          }
          if (matId == -1)
            matId = boundaryCandidate;

          if (matId != -1) {
            multiMesh->insertNextLine(
                {static_cast<unsigned>(boundaryNodes[i].second),
                 static_cast<unsigned>(boundaryNodes[i + 1].second)});
            if (matIds)
              matIds->push_back(matId);
          }
        }
      };

      double minX = totalMinimum[0] * gridDelta;
      double maxX = totalMaximum[0] * gridDelta;
      double minY = totalMinimum[1] * gridDelta;
      double maxY = totalMaximum[1] * gridDelta;

      capBoundary(0, minX, minY, maxY);
      capBoundary(0, maxX, minY, maxY);
      capBoundary(1, minY, minX, maxX);
      capBoundary(1, maxY, minX, maxX);

    } else {
      for (int d = 0; d < D; ++d) {
        for (int side = 0; side < 2; ++side) {
          int boundaryCoord = (side == 0) ? totalMinimum[d] : totalMaximum[d];

          int u = (d + 1) % D;
          int v = (d + 2) % D;
          if (D == 3 && u > v)
            std::swap(u, v);

          int uMin = totalMinimum[u];
          int uMax = totalMaximum[u];
          int vMin = (D == 3) ? totalMinimum[v] : 0;
          int vMax = (D == 3) ? totalMaximum[v] : 0;

          for (int i = uMin; i < uMax; ++i) {
            for (int j = vMin; j < ((D == 3) ? vMax : 1); ++j) {
              double gridCentroid[D];
              gridCentroid[d] = static_cast<double>(boundaryCoord);
              gridCentroid[u] = static_cast<double>(i) + 0.5;
              if (D == 3)
                gridCentroid[v] = static_cast<double>(j) + 0.5;

              int matId = -1;
              int boundaryCandidate = -1;
              viennahrle::Index<D> queryIdx;
              for (int k = 0; k < D; ++k)
                queryIdx[k] = std::round(gridCentroid[k]);

              for (unsigned l = 0; l < levelSets.size(); ++l) {
                viennahrle::ConstDenseIterator<hrleDomainType> it(
                    levelSets[l]->getDomain());
                it.goToIndices(queryIdx);
                double val = it.getValue();
                if (val < -1e-6 * gridDelta) {
                  matId =
                      (materialMap) ? materialMap->getMaterialId(l) : (int)l;
                  break;
                }
                if (val <= 1e-6 * gridDelta && boundaryCandidate == -1) {
                  boundaryCandidate =
                      (materialMap) ? materialMap->getMaterialId(l) : (int)l;
                }
              }
              if (matId == -1)
                matId = boundaryCandidate;

              if (matId != -1) {
                unsigned startNodeId = multiMesh->nodes.size();
                auto matIds = multiMesh->cellData.getScalarData("MaterialIds");

                Vec3D<T> p[4];
                for (int k = 0; k < 4; ++k)
                  p[k][d] = boundaryCoord * gridDelta;

                p[0][u] = i * gridDelta;
                p[0][v] = j * gridDelta;
                p[1][u] = (i + 1) * gridDelta;
                p[1][v] = j * gridDelta;
                p[2][u] = i * gridDelta;
                p[2][v] = (j + 1) * gridDelta;
                p[3][u] = (i + 1) * gridDelta;
                p[3][v] = (j + 1) * gridDelta;

                for (int k = 0; k < 4; ++k)
                  multiMesh->insertNextNode(p[k]);

                multiMesh->insertNextTriangle(
                    {startNodeId, startNodeId + 1, startNodeId + 2});
                multiMesh->insertNextTriangle(
                    {startNodeId + 2, startNodeId + 1, startNodeId + 3});

                if (matIds) {
                  matIds->push_back(matId);
                  matIds->push_back(matId);
                }
              }
            }
          }
        }
      }
    }

    // Convert multiMesh to vtkPolyData
    hullVTK = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (const auto &node : multiMesh->nodes)
      points->InsertNextPoint(node[0], node[1], node[2]);
    hullVTK->SetPoints(points);

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    if constexpr (D == 3) {
      for (const auto &tri : multiMesh->triangles) {
        cells->InsertNextCell(3);
        cells->InsertCellPoint(tri[0]);
        cells->InsertCellPoint(tri[1]);
        cells->InsertCellPoint(tri[2]);
      }
      hullVTK->SetPolys(cells);
    } else {
      for (const auto &line : multiMesh->lines) {
        cells->InsertNextCell(2);
        cells->InsertCellPoint(line[0]);
        cells->InsertCellPoint(line[1]);
      }
      hullVTK->SetLines(cells);
    }

    auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
    if (matIds) {
      vtkSmartPointer<vtkIntArray> materials =
          vtkSmartPointer<vtkIntArray>::New();
      materials->SetName("Material");
      for (auto val : *matIds)
        materials->InsertNextValue(static_cast<int>(val));
      hullVTK->GetCellData()->SetScalars(materials);
    }

    addMetaDataToVTK(hullVTK);

    if (writeToFile) {
      auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName((fileName + "_hull.vtp").c_str());
      writer->SetInputData(hullVTK);
      writer->Write();
    }
  }

public:
  WriteHullMesh() = default;

  WriteHullMesh(SmartPointer<Domain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  /// Level sets wrapping other level sets have to be inserted last.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  /// Set the base file name. "_hull.vtp" will be appended on write.
  void setFileName(std::string passedFileName) {
    fileName = std::move(passedFileName);
  }

  void setWriteToFile(bool passedWriteToFile) {
    writeToFile = passedWriteToFile;
  }

  void setSharpCorners(bool passedSharpCorners) {
    generateSharpCorners = passedSharpCorners;
  }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void setMetaData(const std::unordered_map<std::string, std::vector<double>>
                       &passedMetaData) {
    metaData = passedMetaData;
  }

  void addMetaData(const std::string &key, double value) {
    metaData[key] = std::vector<double>{value};
  }

  void addMetaData(const std::string &key, const std::vector<double> &values) {
    metaData[key] = values;
  }

  void addMetaData(
      const std::unordered_map<std::string, std::vector<double>> &newMetaData) {
    for (const auto &pair : newMetaData)
      metaData[pair.first] = pair.second;
  }

  auto getHullMesh() const { return hullVTK; }

  void apply() { generateSharpHull(); }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(WriteHullMesh)

} // namespace viennals

#endif // VIENNALS_USE_VTK
