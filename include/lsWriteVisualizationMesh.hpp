#ifdef VIENNALS_USE_VTK // this class needs vtk support

#ifndef LS_TO_VISUALIZATION_MESH_HPP
#define LS_TO_VISUALIZATION_MESH_HPP

#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkFloatArray.h>
#include <vtkGeometryFilter.h>
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <hrleDenseIterator.hpp>

#include <lsDomain.hpp>
#include <lsMessage.hpp>

//#define LS_TO_VISUALIZATION_DEBUG
#ifdef LS_TO_VISUALIZATION_DEBUG
#include <vtkXMLRectilinearGridWriter.h>
#endif

/// This algorithm is used to extract tetrahedral volume meshes and triangle
/// hull meshes with material numbers sorted by order of input of level sets. It
/// should ONLY BE USED FOR VISUALIZATION because the algorithm does not
/// guarantee manifold meshes, which should not be a problem for visualization.
/// In order to obtain a hull triangle mesh from the outline of each material,
/// use setExtractHull(true).
template <class T, int D> class lsWriteVisualizationMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;
  using LevelSetsType = std::vector<lsSmartPointer<lsDomain<T, D>>>;
  LevelSetsType levelSets;
  std::string fileName;
  bool extractVolumeMesh = true;
  bool extractHullMesh = false;
  bool bottomRemoved = false;
  static constexpr double LSEpsilon = 1e-2;

  /// This function removes duplicate points and agjusts the pointIDs in the
  /// cells
  /// of a vtkPolyData
  void removeDuplicatePoints(vtkSmartPointer<vtkPolyData> &polyData,
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
  void removeDuplicatePoints(vtkSmartPointer<vtkUnstructuredGrid> &ugrid,
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
  void removeDegenerateTetras(vtkSmartPointer<vtkUnstructuredGrid> &ugrid) {
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
  template <bool removeBottom = false, int gridExtraPoints = 0>
  vtkSmartPointer<vtkRectilinearGrid>
  LS2RectiLinearGrid(lsSmartPointer<lsDomain<T, D>> levelSet,
                     const double LSOffset = 0.,
                     int infiniteMinimum = std::numeric_limits<int>::max(),
                     int infiniteMaximum = -std::numeric_limits<int>::max()) {

    auto &grid = levelSet->getGrid();
    auto &domain = levelSet->getDomain();
    double gridDelta = grid.getGridDelta();
    int numLayers = levelSet->getLevelSetWidth();

    vtkSmartPointer<vtkFloatArray>
        coords[3]; // needs to be 3 because vtk only knows 3D
    int gridMin = 0, gridMax = 0;
    int openJumpDirection = -1; // direction above the open boundary direction

    // fill grid with offset depending on orientation
    for (unsigned i = 0; i < D; ++i) {
      coords[i] = vtkSmartPointer<vtkFloatArray>::New();

      if (grid.getBoundaryConditions(i) ==
          lsDomain<T, D>::BoundaryType::INFINITE_BOUNDARY) {
        // add one to gridMin and gridMax for numerical stability
        gridMin = std::min(domain.getMinRunBreak(i), infiniteMinimum) -
                  1; // choose the smaller number so that for first levelset the
                     // overall minimum can be chosen
        gridMax = std::max(domain.getMaxRunBreak(i), infiniteMaximum) + 1;

        openJumpDirection = i + 1;

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

    // now iterate over grid and fill with LS values
    // initialise iterator over levelset and pointId to start at gridMin
    vtkIdType pointId = 0;
    bool fixBorderPoints = (gridExtraPoints != 0);

    // typename levelSetType::const_iterator_runs it_l(levelSet);
    // use dense iterator to got to every index location
    hrleConstDenseIterator<typename lsDomain<T, D>::DomainType> it(
        levelSet->getDomain());

    // need to save the current position one dimension above open boundary
    // direction, so we can register a jump in the open boundary direction when
    // it occurs, so we can fix the LS value as follows: if remove_bottom==true,
    // we need to flip the sign otherwise the sign stays the same
    int currentOpenIndex =
        it.getIndices()[(openJumpDirection < D) ? openJumpDirection : 0];

    // Make array to store signed distance function
    vtkSmartPointer<vtkFloatArray> signedDistances =
        vtkSmartPointer<vtkFloatArray>::New();
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetName("SignedDistances");

    // iterate until all grid points have a signed distance value
    while ((pointId < rgrid->GetNumberOfPoints())) {
      double p[3];
      rgrid->GetPoint(pointId, p);
      // create index vector
      hrleVectorType<hrleIndexType, D> indices(
          grid.globalCoordinates2GlobalIndices(p));

      // write the corresponding LSValue
      T value;

      // if indices are outside of domain mark point with max value type
      if (grid.isOutsideOfDomain(indices)) {
        fixBorderPoints = true;
        signedDistances->InsertNextValue(
            signedDistances->GetDataTypeValueMax());
      } else {
        // if inside domain just write the correct value
        if (it.getValue() == lsDomain<T, D>::POS_VALUE || it.isFinished()) {
          value = numLayers;
        } else if (it.getValue() == lsDomain<T, D>::NEG_VALUE) {
          value = -numLayers;
        } else {
          value = it.getValue() + LSOffset;
        }

        if (removeBottom) {
          // if we jump from one end of the domain to the other and are not
          // already in the new run, we need to fix the sign of the run
          if (currentOpenIndex != indices[openJumpDirection]) {
            currentOpenIndex = indices[openJumpDirection];
            if (indices >= it.getIndices()) {
              value = -value;
            }
          }
        }

        signedDistances->InsertNextValue(value * gridDelta);
      }

      // move iterator if point was visited
      if (it.isFinished()) { // when iterator is done just fill all the
                             // remaining points
        ++pointId;
      } else {
        while (it.getIndices() <= indices) {
          it.next();
        }

        ++pointId;

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
      }
    }

    // now need to go through again to fix border points, this is done by
    // mapping existing points onto the points outside of the domain according
    // to the correct boundary conditions
    if (fixBorderPoints) {
      pointId = 0;
      while ((pointId < rgrid->GetNumberOfPoints())) {
        if (signedDistances->GetValue(pointId) ==
            signedDistances->GetDataTypeValueMax()) {
          double p[3];
          rgrid->GetPoint(pointId, p);

          // create index vector
          hrleVectorType<hrleIndexType, D> indices(
              grid.globalCoordinates2GlobalIndices(p));

          // vector for mapped point inside domain
          hrleVectorType<hrleIndexType, D> localIndices =
              grid.globalIndices2LocalIndices(indices);

          // now find Id of point we need to take value from
          int originalPointId = 0;
          for (int i = D - 1; i >= 0; --i) {
            originalPointId *=
                coords[i]->GetNumberOfTuples(); // extent in direction
            originalPointId += localIndices[i] - indices[i];
          }
          originalPointId += pointId;

          // now put value of mapped point in global point
          signedDistances->SetValue(pointId,
                                    signedDistances->GetValue(originalPointId));
        }
        ++pointId;
      }
    }

    // Add the SignedDistances to the grid
    rgrid->GetPointData()->SetScalars(signedDistances);

    return rgrid;
  }

public:
  lsWriteVisualizationMesh() {}

  lsWriteVisualizationMesh(lsSmartPointer<lsDomain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  /// Level sets wrapping other level sets have to be inserted last.
  void insertNextLevelSet(lsSmartPointer<lsDomain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  /// Set the name of the file to export. For volume meshes "_volume.vtu" will
  /// be appended, for hull meshes "_hull.vtp".
  void setFileName(std::string passedFileName) { fileName = passedFileName; }

  /// Whether to extract a hull mesh. Defaults to false
  void setExtractHullMesh(bool passedExtractHullMesh) {
    extractHullMesh = passedExtractHullMesh;
  }

  /// Whether to extract a tetra volume mesh. Defaults to true.
  void setExtractVolumeMesh(bool passedExtractVolumeMesh) {
    extractVolumeMesh = passedExtractVolumeMesh;
  }

  void apply() {
    // check if level sets have enough layers
    for (unsigned i = 0; i < levelSets.size(); ++i) {
      if (levelSets[i]->getLevelSetWidth() < 2) {
        lsMessage::getInstance()
            .addWarning(
                "lsWriteVisualizationMesh: Level Set " + std::to_string(i) +
                " should have a width greater than 1! Conversion might fail!")
            .print();
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
            lsDomain<T, D>::BoundaryType::INFINITE_BOUNDARY) {
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
    if (bottomRemoved) {
      topGrid = LS2RectiLinearGrid<true>(levelSets.back(), 0, totalMinimum,
                                         totalMaximum);
    } else {
      topGrid = LS2RectiLinearGrid<false>(levelSets.back(), 0, totalMinimum,
                                          totalMaximum);
    }
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

    materialMeshes.push_back(clipper->GetOutput());
    materialIds.push_back(0);

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

      // create grid of next LS with slight offset and project into current mesh
      vtkSmartPointer<vtkRectilinearGrid> rgrid =
          vtkSmartPointer<vtkRectilinearGrid>::New();
      if (bottomRemoved) {
        rgrid = LS2RectiLinearGrid<true, 1>(*it, -LSEpsilon * counter,
                                            totalMinimum, totalMaximum);
      } else {
        rgrid = LS2RectiLinearGrid<false, 1>(*it, -LSEpsilon * counter,
                                             totalMinimum, totalMaximum);
      }

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

      materialMeshes.rbegin()[0] = insideClipper->GetOutput();
      materialMeshes.push_back(insideClipper->GetClippedOutput());
      materialIds.push_back(counter);

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
        materialNumberArray->InsertNextValue(
            materialIds[materialMeshes.size() - 1 - i]);
      }
      materialMeshes[materialMeshes.size() - 1 - i]->GetCellData()->SetScalars(
          materialNumberArray);

      // delete all point data, so it is not in ouput
      // TODO this includes signed distance information which could be conserved
      // for debugging also includes wheter a cell was vaild for cutting by the
      // grid
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
    auto volumeVTK = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto hullVTK = vtkSmartPointer<vtkPolyData>::New();
    if (extractVolumeMesh) {
      appendFilter->Update();

      // remove degenerate points and remove cells which collapse to zero volume
      // then
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

      // use 1/1000th of grid spacing for contraction of two similar points, so
      // that tetrahedralisation works correctly
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

      auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      writer->SetFileName((fileName + "_volume.vtu").c_str());
      writer->SetInputData(volumeVTK);
      writer->Write();
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

      auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName((fileName + "_hull.vtp").c_str());
      writer->SetInputData(hullVTK);
      writer->Write();
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsWriteVisualizationMesh)

#endif // LS_TO_VISUALIZATION_MESH_HPP
#endif // VIENNALS_USE_VTK