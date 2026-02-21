#pragma once

#ifdef VIENNALS_USE_VTK // this class needs vtk support

#include <vtkAppendFilter.h>
#include <iostream>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkContourTriangulator.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkFloatArray.h>
#include <vtkGeometryFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkClipDataSet.h>
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkStripper.h>
#include <vtkThreshold.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <hrleDenseIterator.hpp>

#include <lsDomain.hpp>
#include <lsMaterialMap.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <unordered_map>
#include <utility>

// #define LS_TO_VISUALIZATION_DEBUG
#ifdef LS_TO_VISUALIZATION_DEBUG
#include <vtkXMLRectilinearGridWriter.h>
#endif

namespace viennals {

using namespace viennacore;

/// This algorithm is used to extract tetrahedral volume meshes and triangle
/// hull meshes with material numbers sorted by order of input of level sets. It
/// should ONLY BE USED FOR VISUALIZATION because the algorithm does not
/// guarantee manifold meshes, which should not be a problem for visualization.
/// In order to obtain a hull triangle mesh from the outline of each material,
/// use setExtractHull(true).
template <class T, int D> class WriteVisualizationMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;
  using LevelSetsType = std::vector<SmartPointer<Domain<T, D>>>;
  LevelSetsType levelSets;
  SmartPointer<MaterialMap> materialMap = nullptr;
  std::string fileName;
  bool extractVolumeMesh = true;
  bool extractHullMesh = false;
  bool bottomRemoved = false;
  bool writeToFile = true;
  bool generateSharpCorners = true;
  double LSEpsilon = 1e-2;
  std::unordered_map<std::string, std::vector<double>> metaData;

  vtkSmartPointer<vtkUnstructuredGrid> volumeVTK;
  vtkSmartPointer<vtkPolyData> hullVTK;

  /// This function removes duplicate points and agjusts the pointIDs in the
  /// cells
  /// of a vtkPolyData
  static void removeDuplicatePoints(vtkSmartPointer<vtkPolyData> &polyData,
                                    const double tolerance) {

    vtkSmartPointer<vtkPolyData> newPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkIncrementalOctreePointLocator> ptInserter =
        vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
    ptInserter->SetTolerance(tolerance);

    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

    // get bounds
    double gridBounds[6];
    polyData->GetBounds(gridBounds);

    // start point insertion from original points
    std::vector<vtkIdType> newPointIds;
    newPointIds.reserve(polyData->GetNumberOfPoints());
    ptInserter->InitPointInsertion(newPoints, gridBounds);

    // make new point list
    for (vtkIdType pointId = 0; pointId < polyData->GetNumberOfPoints();
         ++pointId) {
      vtkIdType globalPtId = 0;
      ptInserter->InsertUniquePoint(polyData->GetPoint(pointId), globalPtId);
      newPointIds.push_back(globalPtId);
    }

    // now set the new points to the unstructured grid
    newPolyData->SetPoints(newPoints);

    // go through all cells and change point ids to match the new ids
    vtkSmartPointer<vtkCellArray> oldCells = polyData->GetPolys();
    vtkSmartPointer<vtkCellArray> newCells =
        vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkIdList> cellPoints = vtkIdList::New();
    oldCells->InitTraversal();
    while (oldCells->GetNextCell(cellPoints)) {
      for (vtkIdType pointId = 0; pointId < cellPoints->GetNumberOfIds();
           ++pointId) {
        cellPoints->SetId(pointId, newPointIds[cellPoints->GetId(pointId)]);
      }
      // insert same cell with new points
      newCells->InsertNextCell(cellPoints);
    }

    newPolyData->SetPolys(newCells);

    // conserve all cell data
    // TODO transfer point data as well (do with "InsertTuples (vtkIdList
    // *dstIds, vtkIdList *srcIds, vtkAbstractArray *source) override)" of
    // vtkDataArray class)
    newPolyData->GetCellData()->ShallowCopy(polyData->GetCellData());

    // set ugrid to the newly created grid
    polyData = newPolyData;
  }

  /// This function removes duplicate points and agjusts the pointIDs in the
  /// cells of a vtkUnstructuredGrid
  static void removeDuplicatePoints(vtkSmartPointer<vtkUnstructuredGrid> &ugrid,
                                    const double tolerance) {

    vtkSmartPointer<vtkUnstructuredGrid> newGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkIncrementalOctreePointLocator> ptInserter =
        vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
    ptInserter->SetTolerance(tolerance);

    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

    // get bounds
    double gridBounds[6];
    ugrid->GetBounds(gridBounds);

    // start point insertion from original points
    std::vector<vtkIdType> newPointIds;
    newPointIds.reserve(ugrid->GetNumberOfPoints());
    ptInserter->InitPointInsertion(newPoints, gridBounds);

    // make new point list
    for (vtkIdType pointId = 0; pointId < ugrid->GetNumberOfPoints();
         ++pointId) {
      vtkIdType globalPtId = 0;
      ptInserter->InsertUniquePoint(ugrid->GetPoint(pointId), globalPtId);
      newPointIds.push_back(globalPtId);
    }

    // now set the new points to the unstructured grid
    newGrid->SetPoints(newPoints);

    // go through all cells and change point ids to match the new ids
    for (vtkIdType cellId = 0; cellId < ugrid->GetNumberOfCells(); ++cellId) {
      vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
      ugrid->GetCellPoints(cellId, cellPoints);
      for (vtkIdType pointId = 0; pointId < cellPoints->GetNumberOfIds();
           ++pointId) {
        cellPoints->SetId(pointId, newPointIds[cellPoints->GetId(pointId)]);
      }
      // insert same cell with new points
      newGrid->InsertNextCell(ugrid->GetCell(cellId)->GetCellType(),
                              cellPoints);
    }

    // conserve all cell data
    // TODO transfer point data as well (do with "InsertTuples (vtkIdList
    // *dstIds, vtkIdList *srcIds, vtkAbstractArray *source) override)" of
    // vtkDataArray class)
    newGrid->GetCellData()->ShallowCopy(ugrid->GetCellData());

    // set ugrid to the newly created grid
    ugrid = newGrid;
  }

  /// This function removes all cells which contain a point more than once
  static void
  removeDegenerateTetras(vtkSmartPointer<vtkUnstructuredGrid> &ugrid) {
    vtkSmartPointer<vtkUnstructuredGrid> newGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    // need to copy material numbers
    vtkSmartPointer<vtkIntArray> materialNumberArray =
        vtkSmartPointer<vtkIntArray>::New();
    materialNumberArray->SetNumberOfComponents(1);
    materialNumberArray->SetName("Material");

    // see if material is defined
    int arrayIndex;
    vtkDataArray *matArray =
        ugrid->GetCellData()->GetArray("Material", arrayIndex);
    const int &materialArrayIndex = arrayIndex;

    // go through all cells and delete those with duplicate entries
    for (vtkIdType cellId = 0; cellId < ugrid->GetNumberOfCells(); ++cellId) {
      vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
      ugrid->GetCellPoints(cellId, cellPoints);
      bool isDuplicate = false;
      for (vtkIdType pointId = 0; pointId < cellPoints->GetNumberOfIds();
           ++pointId) {
        for (vtkIdType nextId = pointId + 1;
             nextId < cellPoints->GetNumberOfIds(); ++nextId) {
          // if they are the same, remove the cell
          if (cellPoints->GetId(pointId) == cellPoints->GetId(nextId))
            isDuplicate = true;
        }
      }
      if (!isDuplicate) {
        // insert same cell if no degenerate points
        newGrid->InsertNextCell(ugrid->GetCell(cellId)->GetCellType(),
                                cellPoints);
        // if material was defined before, use it now
        if (materialArrayIndex >= 0)
          materialNumberArray->InsertNextValue(matArray->GetTuple1(cellId));
      }
    }

    // just take the old points and point data
    newGrid->SetPoints(ugrid->GetPoints());
    newGrid->GetPointData()->ShallowCopy(ugrid->GetPointData());
    // set material cell data
    newGrid->GetCellData()->SetScalars(materialNumberArray);

    ugrid = newGrid;
  }

  // This function takes a levelset and converts it to a vtkRectilinearGrid
  // The full domain contains values, which are capped at numLayers * gridDelta
  // gridExtraPoints layers of grid points are added to the domain according to
  // boundary conditions
  template <int gridExtraPoints = 0>
  vtkSmartPointer<vtkRectilinearGrid>
  LS2RectiLinearGrid(SmartPointer<Domain<T, D>> levelSet, const double LSOffset,
                     int infiniteMinimum = std::numeric_limits<int>::max(),
                     int infiniteMaximum = -std::numeric_limits<int>::max()) {

    auto &grid = levelSet->getGrid();
    auto &domain = levelSet->getDomain();
    double gridDelta = grid.getGridDelta();
    int numLayers = levelSet->getLevelSetWidth();

    vtkSmartPointer<vtkFloatArray>
        coords[3]; // needs to be 3 because vtk only knows 3D
    int gridMin = 0, gridMax = 0;
    // int openJumpDirection = -1; // direction above the open boundary
    // direction

    // fill grid with offset depending on orientation
    for (unsigned i = 0; i < D; ++i) {
      coords[i] = vtkSmartPointer<vtkFloatArray>::New();

      if (grid.getBoundaryConditions(i) ==
          Domain<T, D>::BoundaryType::INFINITE_BOUNDARY) {
        // add one to gridMin and gridMax for numerical stability
        gridMin = std::min(domain.getMinRunBreak(i), infiniteMinimum) -
                  1; // choose the smaller number so that for first levelset the
                     // overall minimum can be chosen
        gridMax = std::max(domain.getMaxRunBreak(i), infiniteMaximum) + 1;

        // openJumpDirection = i + 1;
      } else {
        gridMin = grid.getMinGridPoint(i) - gridExtraPoints;
        gridMax = grid.getMaxGridPoint(i) + gridExtraPoints;
      }

      for (int x = gridMin; x <= gridMax; ++x) {
        coords[i]->InsertNextValue(x * gridDelta);
      }
    }

    // if we work in 2D, just add 1 grid point at origin
    if (D == 2) {
      coords[2] = vtkSmartPointer<vtkFloatArray>::New();
      coords[2]->InsertNextValue(0);
    }

    vtkSmartPointer<vtkRectilinearGrid> rgrid =
        vtkSmartPointer<vtkRectilinearGrid>::New();

    rgrid->SetDimensions(coords[0]->GetNumberOfTuples(),
                         coords[1]->GetNumberOfTuples(),
                         coords[2]->GetNumberOfTuples());
    rgrid->SetXCoordinates(coords[0]);
    rgrid->SetYCoordinates(coords[1]);
    rgrid->SetZCoordinates(coords[2]);

    // bool fixBorderPoints = (gridExtraPoints != 0);

    // need to save the current position one dimension above open boundary
    // direction, so we can register a jump in the open boundary direction when
    // it occurs, so we can fix the LS value as follows: if remove_bottom==true,
    // we need to flip the sign otherwise the sign stays the same
    // int currentOpenIndex =
    //     it.getIndices()[(openJumpDirection < D) ? openJumpDirection : 0];

    // Make array to store signed distance function
    auto const numGridPoints = rgrid->GetNumberOfPoints();
    vtkSmartPointer<vtkFloatArray> signedDistances =
        vtkSmartPointer<vtkFloatArray>::New();
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetNumberOfTuples(numGridPoints);
    signedDistances->SetName("SignedDistances");

#pragma omp parallel
    {
      // use dense iterator to got to every index location
      viennahrle::ConstDenseIterator<typename Domain<T, D>::DomainType> it(
          levelSet->getDomain());

#pragma omp for
      for (vtkIdType pointId = 0; pointId < numGridPoints; ++pointId) {
        // iterate until all grid points have a signed distance value

        double p[3];
        rgrid->GetPoint(pointId, p);
        // create index vector
        viennahrle::Index<D> indices(grid.globalCoordinates2GlobalIndices(p));

        // write the corresponding LSValue
        T value;

        // if indices are outside of domain map to local indices
        if (grid.isOutsideOfDomain(indices)) {
          indices = grid.globalIndices2LocalIndices(indices);
          // fixBorderPoints = true;
          // signedDistances->InsertNextValue(
          // signedDistances->GetDataTypeValueMax());
        }

        it.goToIndices(indices);
        if (it.getValue() == Domain<T, D>::POS_VALUE) {
          value = numLayers;
        } else if (it.getValue() == Domain<T, D>::NEG_VALUE) {
          value = -numLayers;
        } else {
          value = it.getValue() + LSOffset;
        }

        // if (removeBottom) {
        //   // if we jump from one end of the domain to the other and are not
        //   // already in the new run, we need to fix the sign of the run
        //   if (currentOpenIndex != indices[openJumpDirection]) {
        //     currentOpenIndex = indices[openJumpDirection];
        //     if (indices >= it.getIndices()) {
        //       value = -value;
        //     }
        //   }
        // }

        signedDistances->SetValue(pointId, value * gridDelta);
        // signedDistances->InsertNextValue(value * gridDelta);
      }

      // // advance iterator until it is at correct point
      // while (compare(it_l.end_indices(), indices) < 0) {
      //   it_l.next();
      //   if (it_l.is_finished())
      //     break;
      // }
      // // now move iterator with pointId to get to next point
      // switch (compare(it_l.end_indices(), indices)) {
      // case 0:
      //   it_l.next();
      // default:
      //   ++pointId;
      // }
      // }
    }

    // now need to go through again to fix border points, this is done by
    // mapping existing points onto the points outside the domain according
    // to the correct boundary conditions
    // if (fixBorderPoints) {
    //   vtkIdType pointId = 0;
    //   while ((pointId < rgrid->GetNumberOfPoints())) {
    //     if (signedDistances->GetValue(pointId) ==
    //         signedDistances->GetDataTypeValueMax()) {
    //       double p[3];
    //       rgrid->GetPoint(pointId, p);

    //       // create index vector
    //       viennahrle::Index<D> indices(
    //           grid.globalCoordinates2GlobalIndices(p));

    //       // vector for mapped point inside domain
    //       viennahrle::Index<D> localIndices =
    //           grid.globalIndices2LocalIndices(indices);

    //       // now find ID of point we need to take value from
    //       int originalPointId = 0;
    //       for (int i = D - 1; i >= 0; --i) {
    //         originalPointId *=
    //             coords[i]->GetNumberOfTuples(); // extent in direction
    //         originalPointId += localIndices[i] - indices[i];
    //       }
    //       originalPointId += pointId;

    //       // now put value of mapped point in global point
    //       signedDistances->SetValue(pointId,
    //                                 signedDistances->GetValue(originalPointId));
    //     }
    //     ++pointId;
    //   }
    // }

    // Add the SignedDistances to the grid
    rgrid->GetPointData()->SetScalars(signedDistances);

    return rgrid;
  }

  void addMetaDataToVTK(vtkDataSet *data) const {
    if (metaData.empty()) {
      return;
    }

    // add metadata to field data
    vtkSmartPointer<vtkFieldData> fieldData = data->GetFieldData();
    for (const auto &meta : metaData) {
      if (meta.second.empty())
        continue; // skip empty metadata

      vtkSmartPointer<vtkFloatArray> metaDataArray =
          vtkSmartPointer<vtkFloatArray>::New();
      metaDataArray->SetName(meta.first.c_str());
      metaDataArray->SetNumberOfValues(meta.second.size());
      for (size_t i = 0; i < meta.second.size(); ++i) {
        metaDataArray->SetValue(i, meta.second[i]);
      }
      fieldData->AddArray(metaDataArray);
    }
  }

  // // Helper for trilinear/bilinear interpolation
  // T getInterpolatedValue(viennahrle::ConstDenseIterator<hrleDomainType> &it,
  //                        const double *point, double gridDelta) {
  //   int base[D];
  //   double f[D];
  //   for (int i = 0; i < D; ++i) {
  //     double pos = point[i] / gridDelta;
  //     base[i] = static_cast<int>(std::floor(pos));
  //     f[i] = pos - base[i];
  //   }

  //   if constexpr (D == 2) {
  //     auto getVal = [&](int dx, int dy) {
  //       viennahrle::Index<D> idx(base[0] + dx, base[1] + dy);
  //       it.goToIndices(idx);
  //       return it.getValue();
  //     };
  //     return (1 - f[0]) * (1 - f[1]) * getVal(0, 0) +
  //            f[0] * (1 - f[1]) * getVal(1, 0) +
  //            (1 - f[0]) * f[1] * getVal(0, 1) + f[0] * f[1] * getVal(1, 1);
  //   } else {
  //     // Simple nearest neighbor for 3D fallback if needed, 
  //     // but for 2D visualization this function is primarily used.
  //     // Full trilinear can be added if 3D clipping is required.
  //     viennahrle::Index<D> idx;
  //     for(int i=0; i<D; ++i) idx[i] = std::round(point[i]/gridDelta);
  //     it.goToIndices(idx);
  //     return it.getValue();
  //   }
  // }

  void generateSharpMesh() {
    auto &grid = levelSets[0]->getGrid();
    double gridDelta = grid.getGridDelta();

    // Calculate bounds for the domain (handling infinite boundaries)
    int totalMinimum[D];
    int totalMaximum[D];
    for (int i = 0; i < D; ++i) {
      totalMinimum[i] = std::numeric_limits<int>::max();
      totalMaximum[i] = std::numeric_limits<int>::lowest();
    }

    for (auto &ls : levelSets) {
      if (ls->getNumberOfPoints() == 0) continue;
      auto &dom = ls->getDomain();
      auto &grd = ls->getGrid();
      for (int i = 0; i < D; ++i) {
        int minP = (grd.isNegBoundaryInfinite(i))
                       ? dom.getMinRunBreak(i) - 1
                       : grd.getMinBounds(i);
        int maxP = (grd.isPosBoundaryInfinite(i))
                       ? dom.getMaxRunBreak(i) + 1
                       : grd.getMaxBounds(i);
        if (minP < totalMinimum[i])
          totalMinimum[i] = minP;
        if (maxP > totalMaximum[i])
          totalMaximum[i] = maxP;
      }
    }

    // 1. Generate Surface Mesh
    auto multiMesh = SmartPointer<Mesh<T>>::New();
    ToMultiSurfaceMesh<T, D> converter;
    for (auto &ls : levelSets) {
      converter.insertNextLevelSet(ls);
    }
    if (materialMap)
      converter.setMaterialMap(materialMap);
    converter.setMesh(multiMesh);
    converter.setSharpCorners(generateSharpCorners);
    converter.apply();

    // Manually add boundary caps for all domain boundaries
    // This ensures the hull mesh is closed and provides constraints for the
    // volume mesh
    if constexpr (D == 2) {
      auto capBoundary = [&](int boundaryAxis, double fixedValue,
                             double rangeMin, double rangeMax) {
        int varyingAxis = 1 - boundaryAxis;
        std::vector<std::pair<double, unsigned>> boundaryNodes;
        double tolerance = 1e-6 * gridDelta;

        // 1. Find existing nodes on this boundary
        for (unsigned i = 0; i < multiMesh->nodes.size(); ++i) {
          if (std::abs(multiMesh->nodes[i][boundaryAxis] - fixedValue) <
              tolerance) {
            if (multiMesh->nodes[i][varyingAxis] >= rangeMin - tolerance &&
                multiMesh->nodes[i][varyingAxis] <= rangeMax + tolerance) {
              boundaryNodes.push_back({multiMesh->nodes[i][varyingAxis], i});
            }
          }
        }

        // 2. Add corners if not present
        auto addCorner = [&](double val) {
          bool found = false;
          for (const auto &bn : boundaryNodes) {
            if (std::abs(bn.first - val) < tolerance) {
              found = true;
              break;
            }
          }
          if (!found) {
            Vec3D<T> pos;
            pos[boundaryAxis] = fixedValue;
            pos[varyingAxis] = val;
            pos[2] = 0;
            unsigned id = multiMesh->insertNextNode(pos);
            boundaryNodes.push_back({val, id});
          }
        };
        addCorner(rangeMin);
        addCorner(rangeMax);

        // 3. Sort
        std::sort(boundaryNodes.begin(), boundaryNodes.end());

        // 4. Connect
        auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
        for (size_t i = 0; i < boundaryNodes.size() - 1; ++i) {
          // Check midpoint material
          double mid =
              (boundaryNodes[i].first + boundaryNodes[i + 1].first) * 0.5;

          // Check material at centroid
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
        for (int side = 0; side < 2; ++side) { // 0: min, 1: max
          int boundaryCoord = (side == 0) ? totalMinimum[d] : totalMaximum[d];

          // Dimensions to iterate over (the face perpendicular to d)
          int u = (d + 1) % D;
          int v = (d + 2) % D;
          // Ensure u < v for consistent ordering if D=3
          if (D == 3 && u > v)
            std::swap(u, v);

          int uMin = totalMinimum[u];
          int uMax = totalMaximum[u];
          int vMin = (D == 3) ? totalMinimum[v] : 0;
          int vMax = (D == 3) ? totalMaximum[v] : 0;

          for (int i = uMin; i < uMax; ++i) {
            for (int j = vMin; j < ((D == 3) ? vMax : 1); ++j) {
              // Calculate centroid of the boundary element to check material
              double gridCentroid[D];
              gridCentroid[d] = static_cast<double>(boundaryCoord);
              gridCentroid[u] = static_cast<double>(i) + 0.5;
              if (D == 3)
                gridCentroid[v] = static_cast<double>(j) + 0.5;

              // Check material at centroid
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

                // 3D
                // Add 4 nodes for the quad (split into 2 triangles)
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

    // 2. Convert to VTK PolyData (Hull)
    hullVTK = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (const auto &node : multiMesh->nodes) {
      points->InsertNextPoint(node[0], node[1], node[2]);
    }
    hullVTK->SetPoints(points);

    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    if constexpr (D == 3) {
      for (const auto &tri : multiMesh->triangles) {
        polys->InsertNextCell(3);
        polys->InsertCellPoint(tri[0]);
        polys->InsertCellPoint(tri[1]);
        polys->InsertCellPoint(tri[2]);
      }
      hullVTK->SetPolys(polys);
    } else {
      for (const auto &line : multiMesh->lines) {
        polys->InsertNextCell(2);
        polys->InsertCellPoint(line[0]);
        polys->InsertCellPoint(line[1]);
      }
      hullVTK->SetLines(polys);
    }

    // Add Material IDs to Hull
    auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
    if (matIds) {
      vtkSmartPointer<vtkIntArray> materials =
          vtkSmartPointer<vtkIntArray>::New();
      materials->SetName("Material");
      for (auto val : *matIds) {
        materials->InsertNextValue(static_cast<int>(val));
      }
      hullVTK->GetCellData()->SetScalars(materials);
    }

    // Add metadata to hull
    addMetaDataToVTK(hullVTK);

    if (extractHullMesh && writeToFile) {
      auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName((fileName + "_hull.vtp").c_str());
      writer->SetInputData(hullVTK);
      writer->Write();
    }

    // 3. Generate Volume Mesh
    if (extractVolumeMesh) {
      if constexpr (D == 2) {
        vtkSmartPointer<vtkAppendFilter> appendFilter =
            vtkSmartPointer<vtkAppendFilter>::New();

        // Get list of unique material IDs from hullVTK
        vtkDataArray *matIdsArray =
            hullVTK->GetCellData()->GetArray("Material");
        std::set<int> uniqueMatIds;
        if (matIdsArray) {
          for (vtkIdType i = 0; i < matIdsArray->GetNumberOfTuples(); ++i) {
            uniqueMatIds.insert(static_cast<int>(matIdsArray->GetTuple1(i)));
          }
        }

        for (int matId : uniqueMatIds) {
          // Extract lines for this material
          vtkSmartPointer<vtkThreshold> threshold =
              vtkSmartPointer<vtkThreshold>::New();
          threshold->SetInputData(hullVTK);
          threshold->SetInputArrayToProcess(
              0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Material");
          threshold->SetLowerThreshold(matId);
          threshold->SetUpperThreshold(matId);
          threshold->Update();

          // Convert to PolyData for Stripper
          vtkSmartPointer<vtkGeometryFilter> geom =
              vtkSmartPointer<vtkGeometryFilter>::New();
          geom->SetInputData(threshold->GetOutput());
          geom->Update();

          // Create loops
          vtkSmartPointer<vtkStripper> stripper =
              vtkSmartPointer<vtkStripper>::New();
          stripper->SetInputData(geom->GetOutput());
          stripper->JoinContiguousSegmentsOn();
          stripper->Update();

          // Triangulate
          vtkSmartPointer<vtkContourTriangulator> triangulator =
              vtkSmartPointer<vtkContourTriangulator>::New();
          triangulator->SetInputData(stripper->GetOutput());
          triangulator->Update();

          vtkSmartPointer<vtkPolyData> triPoly = triangulator->GetOutput();

          // Set Material ID
          vtkSmartPointer<vtkIntArray> matIdArray =
              vtkSmartPointer<vtkIntArray>::New();
          matIdArray->SetName("Material");
          matIdArray->SetNumberOfComponents(1);
          matIdArray->SetNumberOfTuples(triPoly->GetNumberOfCells());
          for (vtkIdType i = 0; i < triPoly->GetNumberOfCells(); ++i)
            matIdArray->SetValue(i, matId);

          triPoly->GetCellData()->AddArray(matIdArray);

          appendFilter->AddInputData(triPoly);
        }
        appendFilter->Update();
        volumeVTK = appendFilter->GetOutput();
      } else {
      vtkSmartPointer<vtkPoints> internalPoints =
          vtkSmartPointer<vtkPoints>::New();

      for (auto &ls : levelSets) {
        viennahrle::ConstSparseIterator<hrleDomainType> it(ls->getDomain());
        for (; !it.isFinished(); it.next()) {
          if (it.getValue() < 0) {
            auto startIndices = it.getStartIndices();
            for (int k = 0; k < (it.getEndIndices()[0] - startIndices[0]); ++k) {
              auto currentIndices = startIndices;
              currentIndices[0] += k;

              bool onBoundary = false;
              for (int d = 0; d < D; ++d) {
                if (currentIndices[d] == grid.getMinGridPoint(d) ||
                    currentIndices[d] == grid.getMaxGridPoint(d)) {
                  onBoundary = true;
                  break;
                }
              }
              if (it.getValue() < (onBoundary ? -1e-4 * gridDelta : -gridDelta)) {
                Vec3D<T> coord{};
                for (int i = 0; i < D; ++i)
                  coord[i] = currentIndices[i] * gridDelta;
                internalPoints->InsertNextPoint(coord[0], coord[1], coord[2]);
              }
            }
          }
        }
      }

      vtkSmartPointer<vtkUnstructuredGrid> ugrid;

      vtkSmartPointer<vtkPoints> allPoints = vtkSmartPointer<vtkPoints>::New();
      // Add surface points
      for (int i = 0; i < points->GetNumberOfPoints(); ++i) {
        allPoints->InsertNextPoint(points->GetPoint(i));
      }
      // Add internal points
      for (int i = 0; i < internalPoints->GetNumberOfPoints(); ++i) {
        allPoints->InsertNextPoint(internalPoints->GetPoint(i));
      }

      vtkSmartPointer<vtkPolyData> polyPoints =
          vtkSmartPointer<vtkPolyData>::New();
      polyPoints->SetPoints(allPoints);

        // Add 8 corners of the domain to the points
        double minX = totalMinimum[0] * gridDelta;
        double maxX = totalMaximum[0] * gridDelta;
        double minY = totalMinimum[1] * gridDelta;
        double maxY = totalMaximum[1] * gridDelta;
        double minZ = totalMinimum[2] * gridDelta;
        double maxZ = totalMaximum[2] * gridDelta;

        allPoints->InsertNextPoint(minX, minY, minZ);
        allPoints->InsertNextPoint(maxX, minY, minZ);
        allPoints->InsertNextPoint(minX, maxY, minZ);
        allPoints->InsertNextPoint(maxX, maxY, minZ);
        allPoints->InsertNextPoint(minX, minY, maxZ);
        allPoints->InsertNextPoint(maxX, minY, maxZ);
        allPoints->InsertNextPoint(minX, maxY, maxZ);
        allPoints->InsertNextPoint(maxX, maxY, maxZ);

        vtkSmartPointer<vtkDelaunay3D> delaunay =
            vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay->SetInputData(polyPoints);
        delaunay->SetTolerance(0.001 * gridDelta);
        delaunay->Update();
        ugrid = delaunay->GetOutput();

      // Filter Cells
      vtkSmartPointer<vtkUnstructuredGrid> finalVol =
          vtkSmartPointer<vtkUnstructuredGrid>::New();
      finalVol->SetPoints(ugrid->GetPoints());
      finalVol->Allocate(ugrid->GetNumberOfCells());

      vtkSmartPointer<vtkIntArray> volMats =
          vtkSmartPointer<vtkIntArray>::New();
      volMats->SetName("Material");

      const bool useMaterialMap = materialMap != nullptr;

      std::vector<viennahrle::ConstDenseIterator<hrleDomainType>> iterators;
      for (const auto &ls : levelSets) {
        iterators.emplace_back(ls->getDomain());
      }

      std::set<int> presentMaterials;
      for (vtkIdType i = 0; i < ugrid->GetNumberOfCells(); ++i) {
        vtkCell *cell = ugrid->GetCell(i);
        vtkPoints *pts = cell->GetPoints();
        int numPts = pts->GetNumberOfPoints();

        // Determine Material
        int matId = -1;
        for (unsigned l = 0; l < levelSets.size(); ++l) {
          int insideCount = 0;
          auto &grid = levelSets[l]->getGrid();

          for (int j = 0; j < numPts; ++j) {
            double p[3];
            pts->GetPoint(j, p);
            viennahrle::Index<D> idx = grid.globalCoordinates2GlobalIndices(p);
            if (grid.isOutsideOfDomain(idx)) {
              idx = grid.globalIndices2LocalIndices(idx);
            }
            iterators[l].goToIndices(idx);
            // Use a small tolerance to include points slightly outside due to sharp corner reconstruction
            if (iterators[l].getValue() <= 1e-6 * gridDelta) {
              insideCount++;
            }
          }

          if (insideCount > numPts / 2) {
            matId = useMaterialMap ? materialMap->getMaterialId(l) : l;
            break;
          }
        }

        if (matId != -1) {
          finalVol->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
          volMats->InsertNextValue(matId);
          presentMaterials.insert(matId);
        }
      }

      finalVol->GetCellData()->SetScalars(volMats);
      volumeVTK = finalVol;
      }

      // Add metadata
      addMetaDataToVTK(volumeVTK);

      if (extractVolumeMesh && writeToFile) {
        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName((fileName + "_volume.vtu").c_str());
        writer->SetInputData(volumeVTK);
        writer->Write();
      }
    }
  }

public:
  WriteVisualizationMesh() = default;

  WriteVisualizationMesh(SmartPointer<Domain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  /// Level sets wrapping other level sets have to be inserted last.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  /// Set the name of the file to export. For volume meshes "_volume.vtu" will
  /// be appended, for hull meshes "_hull.vtp".
  void setFileName(std::string passedFileName) {
    fileName = std::move(passedFileName);
  }

  /// Whether to extract a hull mesh. Defaults to false
  void setExtractHullMesh(bool passedExtractHullMesh) {
    extractHullMesh = passedExtractHullMesh;
  }

  /// Whether to extract a tetra volume mesh. Defaults to true.
  void setExtractVolumeMesh(bool passedExtractVolumeMesh) {
    extractVolumeMesh = passedExtractVolumeMesh;
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

  void setWrappingLayerEpsilon(double epsilon) { LSEpsilon = epsilon; }

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
    for (const auto &pair : newMetaData) {
      metaData[pair.first] = pair.second;
    }
  }

  auto getVolumeMesh() const { return volumeVTK; }

  auto getHullMesh() const { return hullVTK; }

  void apply() {
    if (generateSharpCorners) {
      generateSharpMesh();
      return;
    }

    // check if level sets have enough layers
    for (unsigned i = 0; i < levelSets.size(); ++i) {
      if (levelSets[i]->getLevelSetWidth() < 2) {
        VIENNACORE_LOG_WARNING(
            "WriteVisualizationMesh: Level Set " + std::to_string(i) +
            " should have a width greater than 1! Conversion might fail!");
      }
    }

    const double gridDelta = levelSets[0]->getGrid().getGridDelta();

    // store volume for each material
    std::vector<vtkSmartPointer<vtkUnstructuredGrid>> materialMeshes;
    std::vector<unsigned> materialIds;

    int totalMinimum = std::numeric_limits<int>::max();
    int totalMaximum = -std::numeric_limits<int>::max();
    for (auto &it : levelSets) {
      if (it->getNumberOfPoints() == 0) {
        continue;
      }
      auto &grid = it->getGrid();
      auto &domain = it->getDomain();
      for (unsigned i = 0; i < D; ++i) {
        if (grid.getBoundaryConditions(i) ==
            Domain<T, D>::BoundaryType::INFINITE_BOUNDARY) {
          totalMinimum = std::min(totalMinimum, domain.getMinRunBreak(i));
          totalMaximum = std::max(totalMaximum, domain.getMaxRunBreak(i));
        }
      }
    }

    // create volume mesh for largest LS
    // Use vtkClipDataSet to slice the grid
    vtkSmartPointer<vtkTableBasedClipDataSet> clipper =
        vtkSmartPointer<vtkTableBasedClipDataSet>::New();
    auto topGrid = vtkSmartPointer<vtkRectilinearGrid>::New();
    // if (bottomRemoved) {
    topGrid =
        LS2RectiLinearGrid(levelSets.back(), 0, totalMinimum, totalMaximum);
    // } else {
    //   topGrid = LS2RectiLinearGrid<false>(levelSets.back(), 0,
    //   totalMinimum,
    //                                       totalMaximum);
    // }
#ifdef LS_TO_VISUALIZATION_DEBUG
    {
      auto gwriter = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
      gwriter->SetFileName("./grid_0.vtr");
      gwriter->SetInputData(topGrid);
      gwriter->Write();
      std::cout << "Wrote grid 0" << std::endl;
    }
#endif
    clipper->SetInputData(topGrid);
    clipper->InsideOutOn();
    clipper->SetValue(0.0);
    clipper->GenerateClippedOutputOn(); // TODO remove
    clipper->Update();

#ifdef LS_TO_VISUALIZATION_DEBUG
    {
      auto gwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      gwriter->SetFileName("./clipped.vtu");
      gwriter->SetInputData(clipper->GetClippedOutput());
      gwriter->Write();
      std::cout << "Wrote clipped" << std::endl;
    }
#endif

    const bool useMaterialMap = materialMap != nullptr;
    materialMeshes.emplace_back(clipper->GetOutput());
    materialIds.push_back(useMaterialMap ? materialMap->getMaterialId(0) : 0);

#ifdef LS_TO_VISUALIZATION_DEBUG
    {
      auto gwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      gwriter->SetFileName("./probed_0.vtu");
      gwriter->SetInputData(materialMeshes.front());
      gwriter->Write();
    }
#endif

    unsigned counter = 1;

    // now cut large volume mesh with all the smaller ones
    for (typename LevelSetsType::const_reverse_iterator it =
             ++levelSets.rbegin();
         it != levelSets.rend(); ++it) {
      if (it->get()->getNumberOfPoints() == 0)
        continue; // ignore empty levelSets

      // create grid of next LS with slight offset and project into current
      // mesh
      vtkSmartPointer<vtkRectilinearGrid> rgrid =
          vtkSmartPointer<vtkRectilinearGrid>::New();
      // if (bottomRemoved) {
      rgrid = LS2RectiLinearGrid<1>(*it, -LSEpsilon * counter, totalMinimum,
                                    totalMaximum);
      // } else {
      //   rgrid = LS2RectiLinearGrid<false, 1>(*it, -LSEpsilon * counter,
      //                                        totalMinimum, totalMaximum);
      // }

#ifdef LS_TO_VISUALIZATION_DEBUG
      {
        vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter =
            vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
        gwriter->SetFileName(
            ("./grid_" + std::to_string(counter) + ".vtr").c_str());
        gwriter->SetInputData(rgrid);
        gwriter->Write();
        std::cout << "Wrote grid " << counter << std::endl;
      }
#endif

      // now transfer implicit values to mesh points
      vtkSmartPointer<vtkProbeFilter> probeFilter =
          vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetInputData(materialMeshes.back()); // last element
      probeFilter->SetSourceData(rgrid);
      probeFilter->Update();

#ifdef LS_TO_VISUALIZATION_DEBUG
      {
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        gwriter->SetFileName(
            ("./probed_" + std::to_string(counter) + ".vtu").c_str());
        gwriter->SetInputData(probeFilter->GetOutput());
        gwriter->Write();
        std::cout << "Wrote unstructured grid " << counter << std::endl;
      }
#endif

      // now clip the mesh and save the clipped as the 1st layer and use the
      // inverse for the next layer clipping Use vtkTabelBasedClipDataSet to
      // slice the grid
      vtkSmartPointer<vtkTableBasedClipDataSet> insideClipper =
          vtkSmartPointer<vtkTableBasedClipDataSet>::New();
      insideClipper->SetInputConnection(probeFilter->GetOutputPort());
      insideClipper->GenerateClippedOutputOn();
      insideClipper->Update();

      materialMeshes.back() = insideClipper->GetOutput();
      materialMeshes.emplace_back(insideClipper->GetClippedOutput());
      auto material = counter;
      if (useMaterialMap)
        material = materialMap->getMaterialId(counter);
      materialIds.push_back(material);

      ++counter;
    }

    vtkSmartPointer<vtkAppendFilter> appendFilter =
        vtkSmartPointer<vtkAppendFilter>::New();

    vtkSmartPointer<vtkAppendPolyData> hullAppendFilter =
        vtkSmartPointer<vtkAppendPolyData>::New();

    for (unsigned i = 0; i < materialMeshes.size(); ++i) {

      // write material number in mesh
      vtkSmartPointer<vtkIntArray> materialNumberArray =
          vtkSmartPointer<vtkIntArray>::New();
      materialNumberArray->SetNumberOfComponents(1);
      materialNumberArray->SetName("Material");
      for (unsigned j = 0;
           j <
           materialMeshes[materialMeshes.size() - 1 - i]->GetNumberOfCells();
           ++j) {
        materialNumberArray->InsertNextValue(materialIds[i]);
      }
      materialMeshes[materialMeshes.size() - 1 - i]->GetCellData()->SetScalars(
          materialNumberArray);

      // delete all point data, so it is not in ouput
      // TODO this includes signed distance information which could be
      // conserved for debugging also includes wheter a cell was vaild for
      // cutting by the grid
      vtkSmartPointer<vtkPointData> pointData =
          materialMeshes[materialMeshes.size() - 1 - i]->GetPointData();
      const int numberOfArrays = pointData->GetNumberOfArrays();
      for (int j = 0; j < numberOfArrays; ++j) {
        pointData->RemoveArray(0); // remove first array until none are left
      }

      // if hull mesh should be exported, create hull for each layer and put
      // them together
      if (extractHullMesh) {
        vtkSmartPointer<vtkGeometryFilter> geoFilter =
            vtkSmartPointer<vtkGeometryFilter>::New();
        geoFilter->SetInputData(materialMeshes[materialMeshes.size() - 1 - i]);
        geoFilter->Update();
        hullAppendFilter->AddInputData(geoFilter->GetOutput());
      }

      appendFilter->AddInputData(materialMeshes[materialMeshes.size() - 1 - i]);
    }

    // do not need tetrahedral volume mesh if we do not print volume
    volumeVTK = vtkSmartPointer<vtkUnstructuredGrid>::New();
    hullVTK = vtkSmartPointer<vtkPolyData>::New();
    if (extractVolumeMesh) {
      appendFilter->Update();

      // remove degenerate points and remove cells which collapse to zero
      // volume then
      volumeVTK = appendFilter->GetOutput();
#ifdef LS_TO_VISUALIZATION_DEBUG
      {
        std::cout << "Before duplicate removal: " << std::endl;
        std::cout << "Points: " << volumeVTK->GetNumberOfPoints() << std::endl;
        std::cout << "Cells: " << volumeVTK->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        gwriter->SetFileName("before_removal.vtu");
        gwriter->SetInputData(appendFilter->GetOutput());
        gwriter->Update();
      }
#endif

      // use 1/1000th of grid spacing for contraction of two similar points,
      // so that tetrahedralisation works correctly
      removeDuplicatePoints(volumeVTK, 1e-3 * gridDelta);

#ifdef LS_TO_VISUALIZATION_DEBUG
      {
        std::cout << "After duplicate removal: " << std::endl;
        std::cout << "Points: " << volumeVTK->GetNumberOfPoints() << std::endl;
        std::cout << "Cells: " << volumeVTK->GetNumberOfCells() << std::endl;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        gwriter->SetFileName("after_removal.vtu");
        gwriter->SetInputData(volumeVTK);
        gwriter->Update();
      }
#endif

      // change all 3D cells into tetras and all 2D cells to triangles
      vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter =
          vtkSmartPointer<vtkDataSetTriangleFilter>::New();
      triangleFilter->SetInputData(volumeVTK);
      triangleFilter->Update();
      volumeVTK = triangleFilter->GetOutput();

      // now that only tetras are left, remove the ones with degenerate points
      removeDegenerateTetras(volumeVTK);

      // add meta data to volume mesh
      addMetaDataToVTK(volumeVTK);

      if (writeToFile) {
        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName((fileName + "_volume.vtu").c_str());
        writer->SetInputData(volumeVTK);
        writer->Write();
      }
    }

    // Now make hull mesh if necessary
    if (extractHullMesh) {
      hullAppendFilter->Update();
      hullVTK = hullAppendFilter->GetOutput();
      // use 1/1000th of grid spacing for contraction of two similar points
      removeDuplicatePoints(hullVTK, 1e-3 * gridDelta);

      vtkSmartPointer<vtkTriangleFilter> hullTriangleFilter =
          vtkSmartPointer<vtkTriangleFilter>::New();
      hullTriangleFilter->SetInputData(hullVTK);
      hullTriangleFilter->Update();

      hullVTK = hullTriangleFilter->GetOutput();

      // add meta data to hull mesh
      addMetaDataToVTK(hullVTK);

      if (writeToFile) {
        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName((fileName + "_hull.vtp").c_str());
        writer->SetInputData(hullVTK);
        writer->Write();
      } else {
        hullTriangleFilter->Update();
      }
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(WriteVisualizationMesh)

} // namespace viennals

#endif // VIENNALS_USE_VTK