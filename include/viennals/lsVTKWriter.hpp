#pragma once

#include <fstream>
#include <string>

#include <lsFileFormats.hpp>
#include <lsMesh.hpp>

#include <utility>
#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

#ifdef VIENNALS_USE_VTK
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif // VIENNALS_USE_VTK

namespace viennals {

using namespace viennacore;

/// Class handling the output of an Mesh<> to VTK file types.
template <class T> class VTKWriter {
  SmartPointer<Mesh<T>> mesh = nullptr;
  FileFormatEnum fileFormat = FileFormatEnum::VTK_AUTO;
  std::string fileName;

#ifdef VIENNALS_USE_VTK
  template <class In, class Out>
  void addDataFromMesh(const In &inData, Out outData) const {
    // now add pointData
    for (unsigned i = 0; i < inData.getScalarDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> pointData =
          vtkSmartPointer<vtkFloatArray>::New();
      pointData->SetNumberOfComponents(1);
      pointData->SetName(inData.getScalarDataLabel(i).c_str());
      auto scalars = *(inData.getScalarData(i));
      for (unsigned j = 0; j < inData.getScalarData(i)->size(); ++j) {
        pointData->InsertNextValue(scalars[j]);
      }
      outData->AddArray(pointData);
    }

    // now add vector data
    for (unsigned i = 0; i < inData.getVectorDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> vectorData =
          vtkSmartPointer<vtkFloatArray>::New();
      vectorData->SetNumberOfComponents(3);
      vectorData->SetName(inData.getVectorDataLabel(i).c_str());
      auto vectors = *(inData.getVectorData(i));
      for (unsigned j = 0; j < inData.getVectorData(i)->size(); ++j) {
        vectorData->InsertNextTuple3(vectors[j][0], vectors[j][1],
                                     vectors[j][2]);
      }
      outData->AddArray(vectorData);
    }
  }
#endif // VIENNALS_USE_VTK

public:
  VTKWriter() = default;

  VTKWriter(SmartPointer<Mesh<T>> passedMesh) : mesh(passedMesh) {}

  VTKWriter(SmartPointer<Mesh<T>> passedMesh, std::string passedFileName)
      : mesh(passedMesh), fileName(std::move(passedFileName)) {}

  VTKWriter(SmartPointer<Mesh<T>> passedMesh, FileFormatEnum passedFormat,
            std::string passedFileName)
      : mesh(passedMesh), fileFormat(passedFormat),
        fileName(std::move(passedFileName)) {}

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  /// set file format for file to write. Defaults to VTK_LEGACY.
  void setFileFormat(FileFormatEnum passedFormat) { fileFormat = passedFormat; }

  /// set file name for file to write
  void setFileName(std::string passedFileName) {
    fileName = std::move(passedFileName);
  }

  void apply() {
    // check mesh
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to VTKWriter. Not writing.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      Logger::getInstance()
          .addWarning("No file name specified for VTKWriter. Not writing.")
          .print();
      return;
    }

    if (fileFormat == FileFormatEnum::VTK_AUTO) {
      auto dotPos = fileName.rfind('.');
      if (dotPos == std::string::npos) {
        fileFormat = FileFormatEnum::VTP;
      } else {
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
                          "passed to VTKWriter. Not writing.")
              .print();
          return;
        }
      }
    }

    // check file format
    switch (fileFormat) {
    case FileFormatEnum::VTK_LEGACY:
      writeVTKLegacy(fileName);
      break;
#ifdef VIENNALS_USE_VTK
    case FileFormatEnum::VTP:
      writeVTP(fileName);
      break;
    case FileFormatEnum::VTU:
      writeVTU(fileName);
      break;
#else
    case FileFormatEnum::VTP:
    case FileFormatEnum::VTU:
      Logger::getInstance()
          .addWarning("VTKWriter was built without VTK support. Falling back "
                      "to VTK_LEGACY.")
          .print();
      writeVTKLegacy(fileName);
      break;
#endif
    default:
      Logger::getInstance()
          .addWarning("No valid file format set for VTKWriter. Not writing.")
          .print();
    }
  }

private:
#ifdef VIENNALS_USE_VTK
  void writeVTP(std::string filename) const {
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to VTKWriter.")
          .print();
      return;
    }

    if (filename.find(".vtp") != filename.size() - 4)
      filename += ".vtp";
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    // Points
    vtkSmartPointer<vtkPoints> polyPoints = vtkSmartPointer<vtkPoints>::New();
    for (auto it = mesh->getNodes().begin(); it != mesh->getNodes().end();
         ++it) {
      polyPoints->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]);
    }
    polyData->SetPoints(polyPoints);

    // Vertices
    if (mesh->vertices.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh->vertices.begin(); it != mesh->vertices.end(); ++it) {
        polyCells->InsertNextCell(1);
        polyCells->InsertCellPoint((*it)[0]);
      }
      polyData->SetVerts(polyCells);
    }

    // Lines
    if (mesh->lines.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh->lines.begin(); it != mesh->lines.end(); ++it) {
        polyCells->InsertNextCell(2);
        for (unsigned i = 0; i < 2; ++i) {
          polyCells->InsertCellPoint((*it)[i]);
        }
      }
      polyData->SetLines(polyCells);
    }

    // Triangles
    if (mesh->triangles.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh->triangles.begin(); it != mesh->triangles.end();
           ++it) {
        polyCells->InsertNextCell(3);
        for (unsigned i = 0; i < 3; ++i) {
          polyCells->InsertCellPoint((*it)[i]);
        }
      }
      polyData->SetPolys(polyCells);
    }

    addDataFromMesh(mesh->pointData, polyData->GetPointData());
    addDataFromMesh(mesh->cellData, polyData->GetCellData());

    vtkSmartPointer<vtkXMLPolyDataWriter> pwriter =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    pwriter->SetFileName(filename.c_str());
    pwriter->SetInputData(polyData);
    pwriter->Write();
  }

  void writeVTU(std::string filename) const {
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to VTKWriter.")
          .print();
      return;
    }

    if (filename.find(".vtu") != filename.size() - 4)
      filename += ".vtu";

    vtkSmartPointer<vtkUnstructuredGrid> uGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto it = mesh->getNodes().begin(); it != mesh->getNodes().end();
         ++it) {
      points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]);
    }
    uGrid->SetPoints(points);

    // Now set all cells
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    std::vector<int> cellTypes;
    cellTypes.reserve(mesh->vertices.size() + mesh->lines.size() +
                      mesh->triangles.size() + mesh->tetras.size() +
                      mesh->hexas.size());

    // Vertices
    if (mesh->vertices.size() > 0) {
      for (auto it = mesh->vertices.begin(); it != mesh->vertices.end(); ++it) {
        cells->InsertNextCell(1);
        cells->InsertCellPoint((*it)[0]);
        cellTypes.push_back(1); // vtk Vertex
      }
    }

    // Lines
    if (mesh->lines.size() > 0) {
      for (auto it = mesh->lines.begin(); it != mesh->lines.end(); ++it) {
        cells->InsertNextCell(2);
        for (unsigned i = 0; i < 2; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(3); // vtk Line
      }
    }

    // Triangles
    if (mesh->triangles.size() > 0) {
      for (auto it = mesh->triangles.begin(); it != mesh->triangles.end();
           ++it) {
        cells->InsertNextCell(3);
        for (unsigned i = 0; i < 3; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(5); // vtk Triangle
      }
    }

    // Tetras
    if (mesh->tetras.size() > 0) {
      for (auto it = mesh->tetras.begin(); it != mesh->tetras.end(); ++it) {
        cells->InsertNextCell(4);
        for (unsigned i = 0; i < 4; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(10); // vtk Tetra
      }
    }

    // Hexas
    if (mesh->hexas.size() > 0) {
      for (auto it = mesh->hexas.begin(); it != mesh->hexas.end(); ++it) {
        cells->InsertNextCell(8);
        for (unsigned i = 0; i < 8; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(12); // vtk Hexahedron
      }
    }

    // set cells
    uGrid->SetCells(&(cellTypes[0]), cells);

    addDataFromMesh(mesh->pointData, uGrid->GetPointData());
    addDataFromMesh(mesh->cellData, uGrid->GetCellData());

    // // now add pointData
    // for (unsigned i = 0; i < mesh->cellData.getScalarDataSize(); ++i) {
    //   vtkSmartPointer<vtkFloatArray> pointData =
    //       vtkSmartPointer<vtkFloatArray>::New();
    //   pointData->SetNumberOfComponents(1);
    //   pointData->SetName(mesh->cellData.getScalarDataLabel(i).c_str());
    //   auto scalars = *(mesh->cellData.getScalarData(i));
    //   for (unsigned j = 0; j < scalars.size(); ++j) {
    //     pointData->InsertNextValue(scalars[j]);
    //   }
    //   uGrid->GetCellData()->AddArray(pointData);
    // }

    // // now add vector data
    // for (unsigned i = 0; i < mesh->cellData.getVectorDataSize(); ++i) {
    //   vtkSmartPointer<vtkFloatArray> vectorData =
    //       vtkSmartPointer<vtkFloatArray>::New();
    //   vectorData->SetNumberOfComponents(3);
    //   vectorData->SetName(mesh->cellData.getVectorDataLabel(i).c_str());
    //   auto vectors = *(mesh->cellData.getVectorData(i));
    //   for (unsigned j = 0; j < vectors.size(); ++j) {
    //     vectorData->InsertNextTuple3(vectors[j][0], vectors[j][1],
    //                                  vectors[j][2]);
    //   }
    //   uGrid->GetCellData()->AddArray(vectorData);
    // }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    owriter->SetFileName(filename.c_str());
    owriter->SetInputData(uGrid);
    owriter->Write();
  }

#endif // VIENNALS_USE_VTK

  void writeVTKLegacy(const std::string &filename) {
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to VTKWriter.")
          .print();
      return;
    }

    std::ofstream f(filename.c_str());

    f << "# vtk DataFile Version 2.0" << std::endl;
    f << ((!mesh->lines.empty()) ? 2 : 3) << "D Surface" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;
    f << "POINTS " << mesh->nodes.size() << " float" << std::endl;

    // print nodes
    for (unsigned int i = 0; i < mesh->nodes.size(); i++) {
      for (int j = 0; j < 3; j++)
        f << static_cast<float>(mesh->nodes[i][j]) << " ";
      f << std::endl;
    }

    const unsigned numberOfCells = mesh->vertices.size() + mesh->lines.size() +
                                   mesh->triangles.size() +
                                   mesh->tetras.size() + mesh->hexas.size();
    const unsigned cellDataSize =
        2 * mesh->vertices.size() + 3 * mesh->lines.size() +
        4 * mesh->triangles.size() + 5 * mesh->tetras.size() +
        9 * mesh->hexas.size();

    f << "CELLS " << numberOfCells << " " << cellDataSize << std::endl;

    // print elements
    for (unsigned int i = 0; i < mesh->vertices.size(); i++) {
      f << 1 << " " << mesh->vertices[i][0] << std::endl;
    }
    for (unsigned int i = 0; i < mesh->lines.size(); i++) {
      f << 2 << " ";
      for (int j = 0; j < 2; j++)
        f << mesh->lines[i][j] << " ";
      f << std::endl;
    }

    for (unsigned int i = 0; i < mesh->triangles.size(); i++) {
      f << 3 << " ";
      for (int j = 0; j < 3; j++)
        f << mesh->triangles[i][j] << " ";
      f << std::endl;
    }

    for (unsigned int i = 0; i < mesh->tetras.size(); i++) {
      f << 4 << " ";
      for (int j = 0; j < 4; j++)
        f << mesh->tetras[i][j] << " ";
      f << std::endl;
    }

    for (unsigned int i = 0; i < mesh->hexas.size(); i++) {
      f << 8 << " ";
      for (int j = 0; j < 8; j++)
        f << mesh->hexas[i][j] << " ";
      f << std::endl;
    }

    f << "CELL_TYPES " << numberOfCells << std::endl;
    for (unsigned i = 0; i < mesh->vertices.size(); ++i)
      f << 1 << std::endl;

    for (unsigned i = 0; i < mesh->lines.size(); ++i)
      f << 3 << std::endl;

    for (unsigned i = 0; i < mesh->triangles.size(); ++i)
      f << 5 << std::endl;

    for (unsigned i = 0; i < mesh->tetras.size(); ++i)
      f << 10 << std::endl;

    for (unsigned i = 0; i < mesh->hexas.size(); ++i)
      f << 12 << std::endl;

    // WRITE POINT DATA
    if (mesh->pointData.getScalarDataSize() ||
        mesh->pointData.getVectorDataSize()) {
      Logger::getInstance()
          .addWarning("Point data output not supported for legacy VTK output. "
                      "Point data is ignored.")
          .print();
    }

    // WRITE SCALAR DATA
    if (mesh->cellData.getScalarDataSize()) {
      f << "CELL_DATA " << mesh->cellData.getScalarData(0)->size() << std::endl;
      for (unsigned i = 0; i < mesh->cellData.getScalarDataSize(); ++i) {
        auto scalars = *(mesh->cellData.getScalarData(i));
        f << "SCALARS " << mesh->cellData.getScalarDataLabel(i) << " float"
          << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for (unsigned j = 0; j < scalars.size(); ++j) {
          f << ((std::abs(scalars[j]) < 1e-6) ? 0.0 : scalars[j]) << std::endl;
        }
      }
    }

    // WRITE VECTOR DATA
    if (mesh->cellData.getVectorDataSize()) {
      if (!mesh->cellData.getScalarDataSize())
        f << "CELL_DATA " << mesh->cellData.getVectorData(0)->size()
          << std::endl;
      for (unsigned i = 0; i < mesh->cellData.getVectorDataSize(); ++i) {
        auto vectors = *(mesh->cellData.getVectorData(i));
        f << "VECTORS " << mesh->cellData.getVectorDataLabel(i) << " float"
          << std::endl;
        for (unsigned j = 0; j < vectors.size(); ++j) {
          for (unsigned k = 0; k < 3; ++k) {
            f << vectors[j][k] << " ";
          }
          f << std::endl;
        }
      }
    }

    f.close();
  }
};

} // namespace viennals
