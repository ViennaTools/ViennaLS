#ifndef LS_VTK_WRITER_HPP
#define LS_VTK_WRITER_HPP

#include <fstream>
#include <string>

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

class lsVTKWriter {
  const lsMesh &mesh;

public:
  lsVTKWriter(lsMesh &passedMesh) : mesh(passedMesh) {}

#ifndef VIENNALS_USE_VTK
  void writeVTP(std::string) const {
    lsMessage::getInstance()
        .addWarning("ViennaLS was built without VTK support. VTK outputs not "
                    "supported.")
        .print();
  }

  void writeVTU(std::string) const { writeVTP(""); }
#else  // VIENNALS_USE_VTK
  void writeVTP(std::string filename) const {
    if (filename.find(".vtp") != filename.size() - 4)
      filename += ".vtp";
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    // Points
    vtkSmartPointer<vtkPoints> polyPoints = vtkSmartPointer<vtkPoints>::New();
    for (auto it = mesh.getNodes().begin(); it != mesh.getNodes().end(); ++it) {
      polyPoints->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]);
    }
    polyData->SetPoints(polyPoints);

    // Vertices
    if (mesh.vertices.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh.vertices.begin(); it != mesh.vertices.end(); ++it) {
        polyCells->InsertNextCell(1);
        polyCells->InsertCellPoint((*it)[0]);
      }
      polyData->SetVerts(polyCells);
    }

    // Lines
    if (mesh.lines.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh.lines.begin(); it != mesh.lines.end(); ++it) {
        polyCells->InsertNextCell(2);
        for (unsigned i = 0; i < 2; ++i) {
          polyCells->InsertCellPoint((*it)[i]);
        }
      }
      polyData->SetLines(polyCells);
    }

    // Triangles
    if (mesh.triangles.size() > 0) {
      vtkSmartPointer<vtkCellArray> polyCells =
          vtkSmartPointer<vtkCellArray>::New();
      for (auto it = mesh.triangles.begin(); it != mesh.triangles.end(); ++it) {
        polyCells->InsertNextCell(3);
        for (unsigned i = 0; i < 3; ++i) {
          polyCells->InsertCellPoint((*it)[i]);
        }
      }
      polyData->SetPolys(polyCells);
    }

    // now add pointData
    for (unsigned i = 0; i < mesh.scalarData.size(); ++i) {
      vtkSmartPointer<vtkFloatArray> pointData =
          vtkSmartPointer<vtkFloatArray>::New();
      pointData->SetNumberOfComponents(1);
      pointData->SetName(mesh.scalarDataLabels[i].c_str());
      for (unsigned j = 0; j < mesh.scalarData[i].size(); ++j) {
        pointData->InsertNextValue(mesh.scalarData[i][j]);
      }
      polyData->GetCellData()->AddArray(pointData);
    }

    // now add vector data
    for (unsigned i = 0; i < mesh.vectorData.size(); ++i) {
      vtkSmartPointer<vtkFloatArray> vectorData =
          vtkSmartPointer<vtkFloatArray>::New();
      vectorData->SetNumberOfComponents(3);
      vectorData->SetName(mesh.vectorDataLabels[i].c_str());
      for (unsigned j = 0; j < mesh.vectorData[i].size(); ++j) {
        vectorData->InsertNextTuple3(mesh.vectorData[i][j][0],
                                     mesh.vectorData[i][j][1],
                                     mesh.vectorData[i][j][2]);
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
    if (filename.find(".vtu") != filename.size() - 4)
      filename += ".vtu";

    vtkSmartPointer<vtkUnstructuredGrid> uGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (auto it = mesh.getNodes().begin(); it != mesh.getNodes().end(); ++it) {
      points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]);
    }
    uGrid->SetPoints(points);

    // Now set all cells
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    std::vector<int> cellTypes;
    cellTypes.reserve(mesh.vertices.size() + mesh.lines.size() +
                      mesh.triangles.size() + mesh.tetras.size() +
                      mesh.hexas.size());

    // Vertices
    if (mesh.vertices.size() > 0) {
      for (auto it = mesh.vertices.begin(); it != mesh.vertices.end(); ++it) {
        cells->InsertNextCell(1);
        cells->InsertCellPoint((*it)[0]);
        cellTypes.push_back(1); // vtk Vertex
      }
    }

    // Lines
    if (mesh.lines.size() > 0) {
      for (auto it = mesh.lines.begin(); it != mesh.lines.end(); ++it) {
        cells->InsertNextCell(2);
        for (unsigned i = 0; i < 2; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(3); // vtk Line
      }
    }

    // Triangles
    if (mesh.triangles.size() > 0) {
      for (auto it = mesh.triangles.begin(); it != mesh.triangles.end(); ++it) {
        cells->InsertNextCell(3);
        for (unsigned i = 0; i < 3; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(5); // vtk Triangle
      }
    }

    // Tetras
    if (mesh.tetras.size() > 0) {
      for (auto it = mesh.tetras.begin(); it != mesh.tetras.end(); ++it) {
        cells->InsertNextCell(4);
        for (unsigned i = 0; i < 4; ++i) {
          cells->InsertCellPoint((*it)[i]);
        }
        cellTypes.push_back(10); // vtk Tetra
      }
    }

    // Hexas
    if (mesh.hexas.size() > 0) {
      for (auto it = mesh.hexas.begin(); it != mesh.hexas.end(); ++it) {
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
    for (unsigned i = 0; i < mesh.scalarData.size(); ++i) {
      vtkSmartPointer<vtkFloatArray> pointData =
          vtkSmartPointer<vtkFloatArray>::New();
      pointData->SetNumberOfComponents(1);
      pointData->SetName(mesh.scalarDataLabels[i].c_str());
      for (unsigned j = 0; j < mesh.scalarData[i].size(); ++j) {
        pointData->InsertNextValue(mesh.scalarData[i][j]);
      }
      uGrid->GetCellData()->AddArray(pointData);
    }

    // now add vector data
    for (unsigned i = 0; i < mesh.vectorData.size(); ++i) {
      vtkSmartPointer<vtkFloatArray> vectorData =
          vtkSmartPointer<vtkFloatArray>::New();
      vectorData->SetNumberOfComponents(3);
      vectorData->SetName(mesh.vectorDataLabels[i].c_str());
      for (unsigned j = 0; j < mesh.vectorData[i].size(); ++j) {
        vectorData->InsertNextTuple3(mesh.vectorData[i][j][0],
                                     mesh.vectorData[i][j][1],
                                     mesh.vectorData[i][j][2]);
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
    std::ofstream f(filename.c_str());

    f << "# vtk DataFile Version 2.0" << std::endl;
    f << ((!mesh.lines.empty()) ? 2 : 3) << "D Surface" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;
    f << "POINTS " << mesh.nodes.size() << " float" << std::endl;

    // print nodes
    for (unsigned int i = 0; i < mesh.nodes.size(); i++) {
      for (int j = 0; j < 3; j++)
        f << static_cast<float>(mesh.nodes[i][j]) << " ";
      f << std::endl;
    }

    // print elements
    if (!mesh.vertices.empty()) {
      f << "CELLS " << mesh.vertices.size() << " " << 2 * mesh.vertices.size()
        << std::endl;
      for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
        f << 1 << " " << mesh.vertices[i][0] << std::endl;
      }
      f << "CELL_TYPES " << mesh.vertices.size() << std::endl;
      for (unsigned i = 0; i < mesh.vertices.size(); ++i)
        f << 1 << std::endl;

    } else if (!mesh.lines.empty()) {
      f << "CELLS " << mesh.lines.size() << " " << 3 * mesh.lines.size()
        << std::endl;
      for (unsigned int i = 0; i < mesh.lines.size(); i++) {
        f << 2 << " ";
        for (int j = 0; j < 2; j++)
          f << mesh.lines[i][j] << " ";
        f << std::endl;
      }
      f << "CELL_TYPES " << mesh.lines.size() << std::endl;
      for (unsigned i = 0; i < mesh.lines.size(); ++i)
        f << 3 << std::endl;

    } else if (!mesh.triangles.empty()) {
      f << "CELLS " << mesh.triangles.size() << " " << 4 * mesh.triangles.size()
        << std::endl;
      for (unsigned int i = 0; i < mesh.triangles.size(); i++) {
        f << 3 << " ";
        for (int j = 0; j < 3; j++)
          f << mesh.triangles[i][j] << " ";
        f << std::endl;
      }
      f << "CELL_TYPES " << mesh.triangles.size() << std::endl;
      for (unsigned i = 0; i < mesh.triangles.size(); ++i)
        f << 5 << std::endl;

    } else if (!mesh.tetras.empty()) {
      f << "CELLS " << mesh.tetras.size() << " " << 5 * mesh.tetras.size()
        << std::endl;
      for (unsigned int i = 0; i < mesh.tetras.size(); i++) {
        f << 4 << " ";
        for (int j = 0; j < 4; j++)
          f << mesh.tetras[i][j] << " ";
        f << std::endl;
      }
      f << "CELL_TYPES " << mesh.tetras.size() << std::endl;
      for (unsigned i = 0; i < mesh.tetras.size(); ++i)
        f << 10 << std::endl;
    } else if (!mesh.hexas.empty()) {
      f << "CELLS " << mesh.hexas.size() << " " << 9 * mesh.hexas.size()
        << std::endl;
      for (unsigned int i = 0; i < mesh.hexas.size(); i++) {
        f << 8 << " ";
        for (int j = 0; j < 8; j++)
          f << mesh.hexas[i][j] << " ";
        f << std::endl;
      }
      f << "CELL_TYPES " << mesh.hexas.size() << std::endl;
      for (unsigned i = 0; i < mesh.hexas.size(); ++i)
        f << 12 << std::endl;
    }

    // WRITE SCALAR DATA
    if (!mesh.scalarData.empty()) {
      f << "CELL_DATA " << mesh.scalarData[0].size() << std::endl;
      for (unsigned i = 0; i < mesh.scalarData.size(); ++i) {
        f << "SCALARS " << mesh.scalarDataLabels[i] << " float" << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for (unsigned j = 0; j < mesh.scalarData[i].size(); ++j) {
          f << mesh.scalarData[i][j] << std::endl;
        }
      }
    }

    // WRITE VECTOR DATA
    if (!mesh.vectorData.empty()) {
      if (mesh.scalarData.empty())
        f << "CELL_DATA " << mesh.vectorData[0].size() << std::endl;
      for (unsigned i = 0; i < mesh.vectorData.size(); ++i) {
        f << "VECTORS " << mesh.vectorDataLabels[i] << " float" << std::endl;
        for (unsigned j = 0; j < mesh.vectorData[i].size(); ++j) {
          for (unsigned k = 0; k < 3; ++k) {
            f << mesh.vectorData[i][j][k] << " ";
          }
          f << std::endl;
        }
      }
    }

    f.close();
  }
};

#endif // LS_VTK_WRITER_HPP
