#pragma once

#include <fstream>
#include <string>
#include <unordered_map>

#include <lsFileFormats.hpp>
#include <lsMesh.hpp>

#include <utility>
#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

#ifdef VIENNALS_USE_VTK
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif // VIENNALS_USE_VTK

namespace viennals {

using namespace viennacore;

/// Class handling the import of VTK file types.
template <class T = double> class VTKReader {
  SmartPointer<Mesh<T>> mesh = nullptr;
  FileFormatEnum fileFormat = FileFormatEnum::VTK_AUTO;
  std::string fileName;
  std::unordered_map<std::string, std::vector<T>> metaData;

  unsigned vtk_nodes_for_cell_type[15] = {0, 1, 0, 2, 0, 3, 0, 0,
                                          4, 4, 4, 8, 8, 6, 5};

  // accepts either vtkPolyData or vtkUnstructuredGrid
  void extractFieldData(vtkDataSet *data) {
    vtkFieldData *fieldData = data->GetFieldData();
    if (!fieldData)
      return;

    for (int i = 0; i < fieldData->GetNumberOfArrays(); ++i) {
      vtkDataArray *array = fieldData->GetArray(i);
      if (!array)
        continue;

      const char *name = array->GetName();
      if (!name)
        continue;

      int numTuples = array->GetNumberOfTuples();
      int numComponents = array->GetNumberOfComponents();

      std::vector<T> values;
      values.reserve(numTuples * numComponents);

      for (int t = 0; t < numTuples; ++t) {
        T tuple[9]; // safe default (up to 9 components)
        array->GetTuple(t, tuple);
        for (int c = 0; c < numComponents; ++c)
          values.push_back(tuple[c]);
      }

      metaData[name] = std::move(values);
    }
  }

public:
  VTKReader() = default;

  VTKReader(SmartPointer<Mesh<T>> passedMesh) : mesh(passedMesh) {}

  VTKReader(SmartPointer<Mesh<T>> passedMesh, std::string passedFileName)
      : mesh(passedMesh), fileName(std::move(passedFileName)) {}

  VTKReader(SmartPointer<Mesh<>> passedMesh, FileFormatEnum passedFormat,
            std::string passedFileName)
      : mesh(passedMesh), fileFormat(passedFormat),
        fileName(std::move(passedFileName)) {}

  /// set the mesh the file should be read into
  void setMesh(SmartPointer<Mesh<>> passedMesh) { mesh = passedMesh; }

  /// set file format for file to read. Defaults to VTK_LEGACY.
  void setFileFormat(FileFormatEnum passedFormat) { fileFormat = passedFormat; }

  /// set file name for file to read
  void setFileName(std::string passedFileName) {
    fileName = std::move(passedFileName);
  }

  auto &getMetaData() { return metaData; }

  void apply() {
    // check mesh
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to VTKReader. Not reading.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      Logger::getInstance()
          .addWarning("No file name specified for VTKReader. Not reading.")
          .print();
      return;
    }

    if (fileFormat == FileFormatEnum::VTK_AUTO) {
      auto dotPos = fileName.rfind('.');
      if (dotPos == std::string::npos) {
        Logger::getInstance()
            .addWarning("No valid file format found based on the file ending "
                        "passed to VTKReader. Not reading.")
            .print();
        return;
      }
      auto ending = fileName.substr(dotPos);
      if (ending == ".vtk") {
        fileFormat = FileFormatEnum::VTK_LEGACY;
      } else if (ending == ".vtp") {
        fileFormat = FileFormatEnum::VTP;
      } else if (ending == ".vtu") {
        fileFormat = FileFormatEnum::VTU;
      } else {
        Logger::getInstance()
            .addWarning("No valid file format found based on the file ending "
                        "passed to VTKReader. Not reading.")
            .print();
        return;
      }
    }

    // check file format
    switch (fileFormat) {
    case FileFormatEnum::VTK_LEGACY:
      readVTKLegacy(fileName);
      break;
#ifdef VIENNALS_USE_VTK
    case FileFormatEnum::VTP:
      readVTP(fileName);
      break;
    case FileFormatEnum::VTU:
      readVTU(fileName);
      break;
#else
    case FileFormatEnum::VTP:
    case FileFormatEnum::VTU:
      Logger::getInstance()
          .addWarning(
              "VTKReader was built without VTK support. Only VTK_LEGACY "
              "can be used. File not read.")
          .print();
#endif
    default:
      Logger::getInstance()
          .addWarning("No valid file format set for VTKReader. Not reading.")
          .print();
    }
  }

private:
#ifdef VIENNALS_USE_VTK
  void readVTP(const std::string &filename) {

    mesh->clear();
    vtkSmartPointer<vtkXMLPolyDataReader> pReader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
    pReader->SetFileName(filename.c_str());
    pReader->Update();

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData = pReader->GetOutput();

    mesh->nodes.resize(polyData->GetNumberOfPoints());
    for (unsigned i = 0; i < mesh->nodes.size(); ++i) {
      std::array<double, 3> coords{};
      polyData->GetPoint(i, coords.data());
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
        std::array<unsigned, 1> cell{};
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
        std::array<unsigned, 2> cell{};
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
        std::array<unsigned, 3> cell{};
        for (unsigned i = 0; i < 3; ++i) {
          cell[i] = pointList->GetId(i);
        }
        mesh->triangles.push_back(cell);
      }
    }

    // get point data
    vtkSmartPointer<vtkPointData> pointData =
        vtkSmartPointer<vtkPointData>::New();
    pointData = polyData->GetPointData();

    for (int i = 0; i < static_cast<int>(pointData->GetNumberOfArrays()); ++i) {
      if (vtkDataArray *dataArray = pointData->GetArray(i);
          dataArray->GetNumberOfComponents() == 1) {
        mesh->pointData.insertNextScalarData(
            typename PointData<T>::ScalarDataType(),
            std::string(pointData->GetArrayName(i)));
        auto &scalars = *(mesh->pointData.getScalarData(i));
        scalars.resize(pointData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          scalars[j] = dataArray->GetTuple1(j);
        }
      } else if (dataArray->GetNumberOfComponents() == 3) {
        mesh->pointData.insertNextVectorData(
            typename PointData<T>::VectorDataType(),
            std::string(pointData->GetArrayName(i)));
        auto &vectors = *(mesh->pointData.getVectorData(i));
        vectors.resize(pointData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          std::array<double, 3> vector{};
          dataArray->GetTuple(j, &(vector[0]));
          vectors[j] = vector;
        }
      }
    }

    // get cell data
    vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
    cellData = polyData->GetCellData();

    for (int i = 0; i < cellData->GetNumberOfArrays(); ++i) {
      vtkDataArray *dataArray = cellData->GetArray(i);
      if (cellData->GetNumberOfComponents() == 1) {
        mesh->cellData.insertNextScalarData(
            typename PointData<T>::ScalarDataType(),
            std::string(cellData->GetArrayName(i)));
        auto &scalars = *(mesh->cellData.getScalarData(i));
        scalars.resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          scalars[j] = dataArray->GetTuple1(j);
        }
      } else if (cellData->GetNumberOfComponents() == 3) {
        mesh->cellData.insertNextVectorData(
            typename PointData<T>::VectorDataType(),
            std::string(cellData->GetArrayName(i)));
        auto &vectors = *(mesh->cellData.getVectorData(i));
        vectors.resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          std::array<double, 3> vector{};
          dataArray->GetTuple(j, &(vector[0]));
          vectors[j] = vector;
        }
      }
    }

    // read meta data
    extractFieldData(polyData);
  }

  void readVTU(const std::string &filename) {

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
      std::array<double, 3> coords{};
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
        std::array<unsigned, 1> vert{};
        vert[0] = pointList->GetId(0);
        mesh->vertices.push_back(vert);
      } break;
      case 3: // line
      {
        std::array<unsigned, 2> elements{};
        for (unsigned j = 0; j < 2; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->lines.push_back(elements);
      } break;
      case 5: // triangle
      {
        std::array<unsigned, 3> elements{};
        for (unsigned j = 0; j < 3; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->triangles.push_back(elements);
      } break;
      case 10: // tetra
      {
        std::array<unsigned, 4> elements{};
        for (unsigned j = 0; j < 4; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->tetras.push_back(elements);
      } break;
      case 12: // hexa
      {
        std::array<unsigned, 8> elements{};
        for (unsigned j = 0; j < 8; ++j) {
          elements[j] = pointList->GetId(j);
        }
        mesh->hexas.push_back(elements);
      } break;
      }
    }

    // get point data
    vtkSmartPointer<vtkPointData> pointData =
        vtkSmartPointer<vtkPointData>::New();
    pointData = ugrid->GetPointData();

    for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
      if (vtkDataArray *dataArray = pointData->GetArray(i);
          dataArray->GetNumberOfComponents() == 1) {
        mesh->pointData.insertNextScalarData(
            typename PointData<T>::ScalarDataType(),
            std::string(pointData->GetArrayName(i)));
        auto &scalars = *(mesh->pointData.getScalarData(i));
        scalars.resize(pointData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          scalars[j] = dataArray->GetTuple1(j);
        }
      } else if (dataArray->GetNumberOfComponents() == 3) {
        mesh->pointData.insertNextVectorData(
            typename PointData<T>::VectorDataType(),
            std::string(pointData->GetArrayName(i)));
        auto &vectors = *(mesh->pointData.getVectorData(i));
        vectors.resize(pointData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          std::array<double, 3> vector{};
          dataArray->GetTuple(j, &(vector[0]));
          vectors[j] = vector;
        }
      }
    }

    // get cell data
    vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
    cellData = ugrid->GetCellData();

    for (int i = 0; i < static_cast<int>(cellData->GetNumberOfArrays()); ++i) {
      vtkDataArray *dataArray = cellData->GetArray(i);
      if (cellData->GetNumberOfComponents() == 1) {
        mesh->cellData.insertNextScalarData(
            typename PointData<T>::ScalarDataType(),
            std::string(cellData->GetArrayName(i)));
        auto &scalars = *(mesh->cellData.getScalarData(i));
        scalars.resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          scalars[j] = dataArray->GetTuple1(j);
        }
      } else if (cellData->GetNumberOfComponents() == 3) {
        mesh->cellData.insertNextVectorData(
            typename PointData<T>::VectorDataType(),
            std::string(cellData->GetArrayName(i)));
        auto &vectors = *(mesh->cellData.getVectorData(i));
        vectors.resize(cellData->GetNumberOfTuples());
        for (unsigned j = 0; j < dataArray->GetNumberOfTuples(); ++j) {
          std::array<double, 3> vector{};
          dataArray->GetTuple(j, &(vector[0]));
          vectors[j] = vector;
        }
      }
    }

    // read meta data
    extractFieldData(ugrid);
  }

#endif // VIENNALS_USE_VTK

  void readVTKLegacy(const std::string &filename) {

    mesh->clear();
    // open geometry file
    std::ifstream f(filename.c_str());
    if (!f)
      Logger::getInstance().addError("Could not open geometry file!");
    std::string temp;

    // Check if geometry is an unstructured grid as required
    while (std::getline(f, temp)) {
      if (temp.find("DATASET") != std::string::npos)
        break;
    }
    if (temp.find("UNSTRUCTURED_GRID") == std::string::npos) {
      Logger::getInstance().addError("DATASET is not an UNSTRUCTURED_GRID!");
    }

    // Find POINTS in file to know number of nodes to read in
    while (std::getline(f, temp)) {
      if (temp.find("POINTS") != std::string::npos)
        break;
    }
    int num_nodes = atoi(&temp[temp.find(' ') + 1]);

    mesh->nodes.resize(num_nodes);

    for (int i = 0; i < num_nodes; i++) {
      double coords[3];

      for (double &coord : coords)
        f >> coord;
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

    int num_elems = atoi(&temp[temp.find(' ') + 1]);

    std::ifstream f_ct(filename.c_str()); // stream to read cell CELL_TYPES
    std::ifstream f_m(
        filename.c_str()); // stream for material numbers if they exist

    // advance to cell types and check if there are the right number
    while (std::getline(f_ct, temp)) {
      if (temp.find("CELL_TYPES") == 0)
        break;
    }
    int num_cell_types = atoi(&temp[temp.find(' ') + 1]);
    // need a cell_type for each cell
    if (num_elems != num_cell_types) {
      Logger::getInstance().addError(
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
    std::vector<VectorType<unsigned int, 4>> elements;
    elements.reserve(num_elems);

    std::vector<double> materials;

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
        case 1: {
          std::array<unsigned, 1> elem{};
          f >> elem[0];
          mesh->template getElements<1>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }
        case 3: {
          std::array<unsigned, 2> elem{};
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->template getElements<2>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }
        case 5: // triangle for 2D
        {
          std::array<unsigned, 3> elem{};
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->template getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 10: // tetra for 3D
        {
          std::array<unsigned, 4> elem{};
          for (unsigned j = 0; j < number_nodes; ++j) {
            f >> elem[j];
          }
          mesh->template getElements<4>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        case 9: // this is a quad, so just plit it into two triangles
        {
          std::array<unsigned, 3> elem{};
          for (unsigned j = 0; j < 3; ++j) {
            f >> elem[j];
          }
          mesh->template getElements<3>().push_back(
              elem); // push the first three nodes as a triangle
          materials.push_back(cell_material);

          f >> elem[1]; // replace middle element to create other triangle
          mesh->template getElements<3>().push_back(elem);
          materials.push_back(cell_material);
          break;
        }

        default:
          std::ostringstream oss;
          oss << "VTK Cell type " << cell_type
              << " is not supported. Cell ignored..." << std::endl;
          Logger::getInstance().addWarning(oss.str()).print();
        }
      } else {
        std::ostringstream oss;
        oss << "INVALID CELL TYPE! Expected number of nodes: " << number_nodes
            << ", Found number of nodes: " << elems_fake
            << "; Ignoring element...";
        Logger::getInstance().addError(oss.str());
        // ignore rest of lines
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }

    mesh->cellData.insertNextScalarData(materials, "Material");

    // Now read Cell Data
    int num_cell_data = 0;
    while (std::getline(f, temp)) {
      if (temp.find("CELL_DATA") != std::string::npos) {
        num_cell_data = atoi(&temp[temp.find(' ') + 1]);
        break;
      }
    }
    std::cout << "Read cell data: " << num_cell_data << std::endl;

    while (!f.eof()) {
      std::cout << "reading scalar data" << std::endl;
      while (std::getline(f, temp)) {
        if (temp.find("SCALARS") != std::string::npos)
          break;
      }
      if (f.eof()) {
        break;
      }

      std::string scalarDataName;
      {
        auto firstS = temp.find(' ') + 1;
        auto secondS = temp.find(' ', firstS + 1);
        scalarDataName = temp.substr(firstS, secondS - firstS);
      }
      std::vector<double> scalarData;

      // consume one line, which defines the lookup table
      std::getline(f, temp);
      if (temp != "LOOKUP_TABLE default") {
        Logger::getInstance()
            .addWarning("Wrong lookup table for VTKLegacy: " + temp)
            .print();
      }

      // now read scalar values
      for (int i = 0; i < num_cell_data; ++i) {
        double data;
        f >> data;
        scalarData.push_back(data);
      }

      mesh->cellData.insertNextScalarData(scalarData, scalarDataName);
    }

    f_ct.close();
    f_m.close();
    f.close();
  }
};

} // namespace viennals
