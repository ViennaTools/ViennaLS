#ifndef LS_VTK_READER_HPP
#define LS_VTK_READER_HPP

#include <fstream>
#include <string>

#include <lsMesh.hpp>
#include <lsMessage.hpp>

#ifdef VIENNALS_USE_VTK
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif // VIENNALS_USE_VTK

/// Class handling the import of VTK file types.
class lsVTKReader {
  lsMesh *mesh = nullptr;

  unsigned vtk_nodes_for_cell_type[15] = {0, 1, 0, 2, 0, 3, 0, 0,
                                          4, 4, 4, 8, 8, 6, 5};

public:
  lsVTKReader(lsMesh &passedMesh) : mesh(&passedMesh) {}

  void setMesh(lsMesh &passedMesh) {
    mesh = &passedMesh;
  }

#ifndef VIENNALS_USE_VTK
  void readVTP(std::string) {
    lsMessage::getInstance()
        .addWarning("ViennaLS was built without VTK support. VTK outputs not "
                    "supported.")
        .print();
  }

  void readVTU(std::string) { readVTP(""); }
#else
  void readVTP(std::string filename) {
    if(mesh == nullptr) {
      lsMessage::getInstance().addWarning("No mesh was passed to lsVTKReader.").print();
      return;
    }

    mesh->clear();
    vtkSmartPointer<vtkXMLPolyDataReader> pReader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
    pReader->SetFileName(filename.c_str());
    pReader->Update();

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData = pReader->GetOutput();

    mesh->nodes.resize(polyData->GetNumberOfPoints());
    for (unsigned i = 0; i < mesh->nodes.size(); ++i) {
      hrleVectorType<double, 3> coords;
      polyData->GetPoint(i, &(coords[0]));
      mesh->nodes[i] = coords;
    }

    vtkSmartPointer<vtkCellArray> cellArray =
        vtkSmartPointer<vtkCellArray>::New();
    // get vertices
    {
      mesh->vertices.reserve(polyData->GetNumberOfVerts());
      cellArray = polyData->GetVerts();
      cellArray->InitTraversal();
      vtkIdList *pointList = vtkIdList::New();
      while (cellArray->GetNextCell(pointList)) {
        hrleVectorType<unsigned, 1> cell;
        cell[0] = pointList->GetId(0);
        mesh->vertices.push_back(cell);
      }
    }

    // get lines
    {
      mesh->lines.reserve(polyData->GetNumberOfLines());
      cellArray = polyData->GetLines();
      cellArray->InitTraversal();
      vtkIdList *pointList = vtkIdList::New();
      while (cellArray->GetNextCell(pointList)) {
        hrleVectorType<unsigned, 2> cell;
        for (unsigned i = 0; i < 2; ++i) {
          cell[i] = pointList->GetId(i);
        }
        mesh->lines.push_back(cell);
      }
    }

    // get triangles
    {
      mesh->triangles.reserve(polyData->GetNumberOfPolys());
      cellArray = polyData->GetPolys();
      cellArray->InitTraversal();
      vtkIdList *pointList = vtkIdList::New();
      while (cellArray->GetNextCell(pointList)) {
        hrleVectorType<unsigned, 3> cell;
        for (unsigned i = 0; i < 3; ++i) {
          cell[i] = pointList->GetId(i);
        }
        mesh->triangles.push_back(cell);
      }
    }

    // get cell data
    vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
    cellData = polyData->GetCellData();

    for (unsigned i = 0;
         i < static_cast<unsigned>(cellData->GetNumberOfArrays()); ++i) {
      vtkDataArray *dataArray;
      dataArray = cellData->GetArray(i);
      if (cellData->GetNumberOfComponents() == 1) {
        mesh->scalarDataLabels.push_back(std::string(cellData->GetArrayName(i)));
        mesh->scalarData.push_back(std::vector<double>());
        mesh->scalarData[i].resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++i) {
          mesh->scalarData[i][j] = dataArray->GetTuple1(j);
        }
      } else if (cellData->GetNumberOfComponents() == 3) {
        mesh->vectorDataLabels.push_back(std::string(cellData->GetArrayName(i)));
        mesh->vectorData.push_back(std::vector<hrleVectorType<double, 3>>());
        mesh->vectorData[i].resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++i) {
          hrleVectorType<double, 3> vector;
          dataArray->GetTuple(j, &(vector[0]));
          mesh->vectorData[i][j] = vector;
        }
      }
    }
  }

  void readVTU(std::string filename) {
    if(mesh == nullptr) {
      lsMessage::getInstance().addWarning("No mesh was passed to lsVTKReader.").print();
      return;
    }

    mesh->clear();

    vtkSmartPointer<vtkXMLUnstructuredGridReader> greader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    greader->SetFileName(filename.c_str());
    greader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid = greader->GetOutput();

    // get all points
    mesh->nodes.resize(ugrid->GetNumberOfPoints());
    for (unsigned i = 0; i < mesh->nodes.size(); ++i) {
      hrleVectorType<double, 3> coords;
      ugrid->GetPoint(i, &(coords[0]));
      mesh->nodes[i] = coords;
    }

    // get cells
    for (unsigned i = 0; i < ugrid->GetNumberOfCells(); ++i) {
      vtkIdList *pointList = vtkIdList::New();
      ugrid->GetCellPoints(i, pointList);

      switch (ugrid->GetCellType(i)) {
      case 1: // vert
      {
        hrleVectorType<unsigned, 1> vert;
        vert[0] = pointList->GetId(0);
        mesh->vertices.push_back(vert);
      } break;
      case 3: // line
      {
        hrleVectorType<unsigned, 2> elements;
        for (unsigned j = 0; j < 2; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->lines.push_back(elements);
      } break;
      case 5: // triangle
      {
        hrleVectorType<unsigned, 3> elements;
        for (unsigned j = 0; j < 3; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->triangles.push_back(elements);
      } break;
      case 10: // tetra
      {
        hrleVectorType<unsigned, 4> elements;
        for (unsigned j = 0; j < 4; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->tetras.push_back(elements);
      } break;
      case 12: // hexa
      {
        hrleVectorType<unsigned, 8> elements;
        for (unsigned j = 0; j < 8; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->hexas.push_back(elements);
      } break;
      }
    }

    // get cell data
    vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
    cellData = ugrid->GetCellData();

    for (unsigned i = 0;
         i < static_cast<unsigned>(cellData->GetNumberOfArrays()); ++i) {
      vtkDataArray *dataArray;
      dataArray = cellData->GetArray(i);
      if (cellData->GetNumberOfComponents() == 1) {
        mesh->scalarDataLabels.push_back(std::string(cellData->GetArrayName(i)));
        mesh->scalarData.push_back(std::vector<double>());
        mesh->scalarData[i].resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++i) {
          mesh->scalarData[i][j] = dataArray->GetTuple1(j);
        }
      } else if (cellData->GetNumberOfComponents() == 3) {
        mesh->vectorDataLabels.push_back(std::string(cellData->GetArrayName(i)));
        mesh->vectorData.push_back(std::vector<hrleVectorType<double, 3>>());
        mesh->vectorData[i].resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++i) {
          hrleVectorType<double, 3> vector;
          dataArray->GetTuple(j, &(vector[0]));
          mesh->vectorData[i][j] = vector;
        }
      }
    }
  }
#endif

  void readVTKLegacy(std::string filename) {
    if(mesh == nullptr) {
      lsMessage::getInstance().addWarning("No mesh was passed to lsVTKReader.").print();
      return;
    }

    mesh->clear();
    // open geometry file
    std::ifstream f(filename.c_str());
    if (!f)
      lsMessage::getInstance().addError("Could not open geometry file!");
    std::string temp;

    // Check if geometry is an unstructured grid as required
    while (std::getline(f, temp)) {
      if (temp.find("DATASET") != std::string::npos)
        break;
    }
    if (temp.find("UNSTRUCTURED_GRID") == std::string::npos) {
      lsMessage::getInstance().addError("DATASET is not an UNSTRUCTURED_GRID!");
    }

    // Find POINTS in file to know number of nodes to read in
    while (std::getline(f, temp)) {
      if (temp.find("POINTS") != std::string::npos)
        break;
    }
    int num_nodes = atoi(&temp[temp.find(" ") + 1]);

    mesh->nodes.resize(num_nodes);

    for (int i = 0; i < num_nodes; i++) {
      double coords[3];

      for (int j = 0; j < 3; j++)
        f >> coords[j];
      for (int j = 0; j < 3; j++) {
        mesh->nodes[i][j] = coords[j];
        // int shift_size = shift.size();
        // if (shift_size > j)
        //   mesh->nodes[i][j] += shift[j]; // Assign desired shift
        // if (InputTransformationSigns[j])
        //   mesh->nodes[i][j] = -mesh->nodes[i][j]; // Assign sign
        //   transformation, if needed
        // mesh->nodes[i][j] *= scale; // Scale the geometry according to
        // parameters file
      }
    }

    while (std::getline(f, temp)) {
      if (temp.find("CELLS") == 0)
        break;
    }

    int num_elems = atoi(&temp[temp.find(" ") + 1]);

    std::ifstream f_ct(filename.c_str()); // stream to read cell CELL_TYPES
    std::ifstream f_m(
        filename.c_str()); // stream for material numbers if they exist

    // advance to cell types and check if there are the right number
    while (std::getline(f_ct, temp)) {
      if (temp.find("CELL_TYPES") == 0)
        break;
    }
    int num_cell_types = atoi(&temp[temp.find(" ") + 1]);
    // need a cell_type for each cell
    if (num_elems != num_cell_types) {
      lsMessage::getInstance().addError(
          "Corrupt input geometry! Number of CELLS and CELL_TYPES "
          "is different!");
    }

    bool is_material = true;
    // advance to material if it is specified
    while (std::getline(f_m, temp)) {
      if (temp.find("CELL_DATA") != std::string::npos) {
        std::getline(f_m, temp);
        if ((temp.find("SCALARS material") != std::string::npos) ||
            (temp.find("SCALARS Material") != std::string::npos)) {
          std::getline(f_m, temp);
          break;
        }
      }
    }
    if (f_m.eof()) {
      is_material = false;
    }

    // maximum cell to read in is a tetra with 4 points
    std::vector<hrleVectorType<unsigned int, 4>> elements;
    elements.reserve(num_elems);

    auto materials = mesh->getMaterials();
    materials.clear();

    unsigned elems_fake;
    unsigned cell_type;
    unsigned cell_material;
    for (int i = 0; i < num_elems; i++) {
      f >> elems_fake;
      f_ct >> cell_type;
      if (is_material)
        f_m >> cell_material;
      else
        cell_material =
            1; // if there are no materials specified make all the same

      // check if the correct number of nodes for cell_type is given
      unsigned number_nodes = vtk_nodes_for_cell_type[cell_type];
      if (number_nodes == elems_fake || number_nodes == 0) {
        // check for different types to subdivide them into supported types
        switch (cell_type) {
        case 3: {
          hrleVectorType<unsigned, 2> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->getElements<2>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }
        case 5: // triangle for 2D
        {
          hrleVectorType<unsigned, 3> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 10: // tetra for 3D
        {
          hrleVectorType<unsigned, 4> elem;
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->getElements<4>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 9: // this is a quad, so just plit it into two triangles
        {
          hrleVectorType<unsigned, 3> elem;
          for (unsigned j = 0; j < 3; ++j) {
            f >> elem[j];
          }
          mesh->getElements<3>().push_back(
              elem); // push the first three nodes as a triangle
          materials.push_back(cell_material);

          f >> elem[1]; // replace middle element to create other triangle
          mesh->getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        default:
          std::ostringstream oss;
          oss << "VTK Cell type " << cell_type
              << " is not supported. Cell ignored..." << std::endl;
          lsMessage::getInstance().addWarning(oss.str()).print();
        }
      } else {
        std::ostringstream oss;
        oss << "INVALID CELL TYPE! Expected number of nodes: " << number_nodes
            << ", Found number of nodes: " << elems_fake
            << "; Ignoring element...";
        lsMessage::getInstance().addError(oss.str());
        // ignore rest of lines
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }

    f_ct.close();
    f_m.close();
    f.close();
  }
};

#endif // LS_VTK_READER_HPP
