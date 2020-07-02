#ifndef LS_VTK_WRITER_HPP
#define LS_VTK_WRITER_HPP

#include <fstream>
#include <string>

#include <lsFileFormats.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

#ifdef VIENNALS_USE_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif // VIENNALS_USE_VTK

/// Class handling the output of an lsMesh to VTK file types.
class lsVTKWriter {
  lsSmartPointer<lsMesh> mesh = nullptr;
  lsFileFormatEnum fileFormat = lsFileFormatEnum::VTK_LEGACY;
  std::string fileName;

public:
  lsVTKWriter() {}

  lsVTKWriter(lsSmartPointer<lsMesh> passedMesh) : mesh(passedMesh) {}

  lsVTKWriter(lsSmartPointer<lsMesh> passedMesh, std::string passedFileName)
      : mesh(passedMesh), fileName(passedFileName) {}

  lsVTKWriter(lsSmartPointer<lsMesh> passedMesh, lsFileFormatEnum passedFormat,
              std::string passedFileName)
      : mesh(passedMesh), fileFormat(passedFormat), fileName(passedFileName) {}

  void setMesh(lsSmartPointer<lsMesh> passedMesh) { mesh = passedMesh; }

  /// set file format for file to write. Defaults to VTK_LEGACY.
  void setFileFormat(lsFileFormatEnum passedFormat) {
    fileFormat = passedFormat;
  }

  /// set file name for file to write
  void setFileName(std::string passedFileName) { fileName = passedFileName; }

  void apply() {
    // check mesh
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsVTKWriter. Not writing.")
          .print();
      return;
    }
    // check filename
    if (fileName.empty()) {
      lsMessage::getInstance()
          .addWarning("No file name specified for lsVTKWriter. Not writing.")
          .print();
      return;
    }

    // check file format
    switch (fileFormat) {
    case lsFileFormatEnum::VTK_LEGACY:
      writeVTKLegacy(fileName);
      break;
#ifdef VIENNALS_USE_VTK
    case lsFileFormatEnum::VTP:
      writeVTP(fileName);
      break;
    case lsFileFormatEnum::VTU:
      writeVTU(fileName);
      break;
#else
    case lsFileFormatEnum::VTP:
    case lsFileFormatEnum::VTU:
      lsMessage::getInstance()
          .addWarning(
              "lsVTKWriter was built without VTK support. Only VTK_LEGACY "
              "can be used. File not written.")
          .print();
#endif
    default:
      lsMessage::getInstance()
          .addWarning("No valid file format set for lsVTKWriter. Not writing.")
          .print();
    }
  }

private:
#ifdef VIENNALS_USE_VTK
  void writeVTP(std::string filename) const {
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsVTKWriter.")
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

    // now add pointData
    for (unsigned i = 0; i < mesh->getScalarDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> pointData =
          vtkSmartPointer<vtkFloatArray>::New();
      pointData->SetNumberOfComponents(1);
      pointData->SetName(mesh->getScalarDataLabel(i).c_str());
      auto scalars = *(mesh->getScalarData(i));
      for (unsigned j = 0; j < mesh->getScalarData(i)->size(); ++j) {
        pointData->InsertNextValue(scalars[j]);
      }
      polyData->GetCellData()->AddArray(pointData);
    }

    // now add vector data
    for (unsigned i = 0; i < mesh->getVectorDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> vectorData =
          vtkSmartPointer<vtkFloatArray>::New();
      vectorData->SetNumberOfComponents(3);
      vectorData->SetName(mesh->getVectorDataLabel(i).c_str());
      auto vectors = *(mesh->getVectorData(i));
      for (unsigned j = 0; j < mesh->getVectorData(i)->size(); ++j) {
        vectorData->InsertNextTuple3(vectors[j][0], vectors[j][1],
                                     vectors[j][2]);
      }
      polyData->GetCellData()->AddArray(vectorData);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> pwriter =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    pwriter->SetFileName(filename.c_str());
    pwriter->SetInputData(polyData);
    pwriter->Write();
  }

  void writeVTU(std::string filename) const {
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsVTKWriter.")
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

    // now add pointData
    for (unsigned i = 0; i < mesh->getScalarDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> pointData =
          vtkSmartPointer<vtkFloatArray>::New();
      pointData->SetNumberOfComponents(1);
      pointData->SetName(mesh->getScalarDataLabel(i).c_str());
      auto scalars = *(mesh->getScalarData(i));
      for (unsigned j = 0; j < scalars.size(); ++j) {
        pointData->InsertNextValue(scalars[j]);
      }
      uGrid->GetCellData()->AddArray(pointData);
    }

    // now add vector data
    for (unsigned i = 0; i < mesh->getVectorDataSize(); ++i) {
      vtkSmartPointer<vtkFloatArray> vectorData =
          vtkSmartPointer<vtkFloatArray>::New();
      vectorData->SetNumberOfComponents(3);
      vectorData->SetName(mesh->getVectorDataLabel(i).c_str());
      auto vectors = *(mesh->getVectorData(i));
      for (unsigned j = 0; j < vectors.size(); ++j) {
        vectorData->InsertNextTuple3(vectors[j][0], vectors[j][1],
                                     vectors[j][2]);
      }
      uGrid->GetCellData()->AddArray(vectorData);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    owriter->SetFileName(filename.c_str());
    owriter->SetInputData(uGrid);
    owriter->Write();
  }

#endif // VIENNALS_USE_VTK

  void writeVTKLegacy(std::string filename) {
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsVTKWriter.")
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

    // WRITE SCALAR DATA
    if (mesh->getScalarDataSize()) {
      f << "CELL_DATA " << mesh->getScalarData(0)->size() << std::endl;
      for (unsigned i = 0; i < mesh->getScalarDataSize(); ++i) {
        auto scalars = *(mesh->getScalarData(i));
        f << "SCALARS " << mesh->getScalarDataLabel(i) << " float" << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for (unsigned j = 0; j < scalars.size(); ++j) {
          f << scalars[j] << std::endl;
        }
      }
    }

    // WRITE VECTOR DATA
    if (mesh->getVectorDataSize()) {
      if (!mesh->getScalarDataSize())
        f << "CELL_DATA " << mesh->getVectorData(0)->size() << std::endl;
      for (unsigned i = 0; i < mesh->getVectorDataSize(); ++i) {
        auto vectors = *(mesh->getVectorData(i));
        f << "VECTORS " << mesh->getVectorDataLabel(i) << " float" << std::endl;
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

#endif // LS_VTK_WRITER_HPP
