#pragma once

#ifdef VIENNALS_VTK_RENDERING

#include <vtkActor.h>
#include <vtkAutoInit.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingUI);

namespace viennals {

template <typename T> class VTKRenderWindow {
public:
  VTKRenderWindow() {
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(800, 600);
    renderWindow->SetWindowName("ViennaLS Render Window");
    renderWindow->SetPosition(100, 100);
    renderWindow->AddRenderer(renderer);
    renderer->SetBackground(0.1, 0.2, 0.3); // Gray background
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);
  }

  ~VTKRenderWindow() {
    if (interactor) {
      interactor->SetRenderWindow(nullptr);
    }
    if (renderWindow) {
      renderWindow->RemoveRenderer(renderer);
    }
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) {
    mesh = passedMesh;
    updatePolyData();
  }

  void render() {
    if (mesh == nullptr) {
      std::cout << "No mesh set for rendering." << std::endl;
      return;
    }
    std::cout << "Rendering mesh with " << mesh->getNodes().size() << " nodes, "
              << mesh->lines.size() << " lines, " << mesh->triangles.size()
              << " triangles." << std::endl;
    renderWindow->Render();
    interactor->Initialize();
    interactor->Start();
  }

  vtkSmartPointer<vtkRenderWindow> getRenderWindow() { return renderWindow; }

private:
  SmartPointer<Mesh<T>> mesh = nullptr;
  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkRenderWindowInteractor> interactor;
  vtkSmartPointer<vtkPolyData> polyData;
  vtkSmartPointer<vtkPolyDataMapper> mapper;
  vtkSmartPointer<vtkActor> actor;

  void updatePolyData() {
    if (mesh == nullptr) {
      return;
    }

    // Remove previous actor if exists
    if (actor) {
      renderer->RemoveActor(actor);
    }

    polyData = vtkSmartPointer<vtkPolyData>::New();

    // Points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (const auto &node : mesh->getNodes()) {
      points->InsertNextPoint(node[0], node[1], node[2]);
    }
    polyData->SetPoints(points);

    // Debug: print bounds
    double bounds[6];
    polyData->GetBounds(bounds);
    std::cout << "PolyData bounds: " << bounds[0] << " " << bounds[1] << " "
              << bounds[2] << " " << bounds[3] << " " << bounds[4] << " "
              << bounds[5] << std::endl;

    // Vertices
    if (!mesh->vertices.empty()) {
      vtkSmartPointer<vtkCellArray> verts =
          vtkSmartPointer<vtkCellArray>::New();
      for (const auto &vertex : mesh->vertices) {
        verts->InsertNextCell(1);
        verts->InsertCellPoint(vertex[0]);
      }
      polyData->SetVerts(verts);
    }

    // Lines
    if (!mesh->lines.empty()) {
      std::cout << "Adding " << mesh->lines.size() << " lines to VTK polydata."
                << std::endl;
      vtkSmartPointer<vtkCellArray> lines =
          vtkSmartPointer<vtkCellArray>::New();
      for (const auto &line : mesh->lines) {
        lines->InsertNextCell(2);
        lines->InsertCellPoint(line[0]);
        lines->InsertCellPoint(line[1]);
      }
      polyData->SetLines(lines);
    }

    // Triangles
    if (!mesh->triangles.empty()) {
      std::cout << "Adding " << mesh->triangles.size()
                << " triangles to VTK polydata." << std::endl;
      vtkSmartPointer<vtkCellArray> polys =
          vtkSmartPointer<vtkCellArray>::New();
      for (const auto &triangle : mesh->triangles) {
        polys->InsertNextCell(3);
        polys->InsertCellPoint(triangle[0]);
        polys->InsertCellPoint(triangle[1]);
        polys->InsertCellPoint(triangle[2]);
      }
      polyData->SetPolys(polys);
    }

    std::cout << "Number of cells in polyData: " << polyData->GetNumberOfCells()
              << std::endl;

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 0.0, 0.0); // Red color
    actor->GetProperty()->SetLineWidth(3.0);       // Thicker lines

    renderer->AddActor(actor);
    renderer->ResetCamera();
  }
};

} // namespace viennals

#endif // VIENNALS_VTK_RENDERING