#pragma once

#ifdef VIENNALS_VTK_RENDERING

#include <lsMaterialMap.hpp>
#include <lsMesh.hpp>

#include <vtkActor.h>
#include <vtkAutoInit.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnstructuredGrid.h>

#ifndef VIENNALS_VTK_MODULE_INIT
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkRenderingUI);
#endif

class ImagePanInteractorStyle : public vtkInteractorStyleImage {
public:
  static ImagePanInteractorStyle *New();
  vtkTypeMacro(ImagePanInteractorStyle, vtkInteractorStyleImage);

  void OnLeftButtonDown() override { this->StartPan(); }

  void OnLeftButtonUp() override { this->EndPan(); }
};

vtkStandardNewMacro(ImagePanInteractorStyle);

namespace viennals {

template <typename T> class VTKRenderWindow {
public:
  VTKRenderWindow() { initialize(); }
  VTKRenderWindow(SmartPointer<Mesh<T>> passedMesh) {
    initialize();
    setMesh(passedMesh);
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
    auto matIds = mesh->getCellData().getScalarData("MaterialIds", false);
    if (matIds) {
      materialIds = *matIds;
    }
    updatePolyData();
  }

  void setVolumeMesh(vtkSmartPointer<vtkUnstructuredGrid> volumeVTK) {
    volumeMesh = volumeVTK;
    updateVolumeMesh();
  }

  void setMaterialIds(const std::vector<T> &ids) { materialIds = ids; }

  auto setBackgroundColor(const std::array<double, 3> &color) {
    backgroundColor = color;
    if (renderer) {
      renderer->SetBackground(backgroundColor.data());
    }

    return *this;
  }

  void render() {
    if (mesh == nullptr && volumeMesh == nullptr) {
      VIENNACORE_LOG_WARNING("No mesh set for rendering.");
      return;
    }

    // Re-attach renderer to the style right before starting, in case the style
    // was changed in the meantime
    if (auto style = interactor->GetInteractorStyle()) {
      style->SetDefaultRenderer(renderer);
      style->SetCurrentRenderer(renderer);
    }
    interactor->SetRenderWindow(renderWindow);
    renderWindow->AddRenderer(renderer);

    renderWindow->Render();
    interactor->Initialize();
    interactor->Start();
  }

  vtkSmartPointer<vtkRenderWindow> getRenderWindow() { return renderWindow; }

  auto enable2DMode() {
    assert(renderer && "Renderer not initialized");
    vtkCamera *cam = renderer->GetActiveCamera();
    cam->ParallelProjectionOn();

    auto style = vtkSmartPointer<ImagePanInteractorStyle>::New();
    // Make sure both the interactor and the style know which renderer to use
    style->SetDefaultRenderer(renderer);
    style->SetCurrentRenderer(renderer);
    interactor->SetInteractorStyle(style);
    return *this;
  }

private:
  void initialize() {
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(backgroundColor.data());

    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetWindowName("ViennaLS Render Window");
    renderWindow->SetSize(windowSize.data());

    // Initialize interactor
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    // Set interactor style
    auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    interactor->SetInteractorStyle(style);

    renderWindow->SetInteractor(interactor);
  }

  void updatePolyData() {
    if (mesh == nullptr) {
      return;
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
      vtkSmartPointer<vtkCellArray> lines =
          vtkSmartPointer<vtkCellArray>::New();
      for (const auto &line : mesh->lines) {
        lines->InsertNextCell(2);
        lines->InsertCellPoint(line[0]);
        lines->InsertCellPoint(line[1]);
      }
      polyData->SetLines(lines);

      enable2DMode();
    }

    // Triangles
    if (!mesh->triangles.empty()) {
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

    // Material IDs as cell data
    bool useMaterialIds =
        !materialIds.empty() &&
        (materialIds.size() == mesh->lines.size() + mesh->triangles.size());
    int minId = std::numeric_limits<int>::max();
    int maxId = std::numeric_limits<int>::min();
    if (useMaterialIds) {
      vtkSmartPointer<vtkIntArray> matIdArray =
          vtkSmartPointer<vtkIntArray>::New();
      matIdArray->SetName("MaterialIds");
      for (const auto &id : materialIds) {
        int mId = static_cast<int>(id);
        matIdArray->InsertNextValue(mId);
        minId = std::min(minId, mId);
        maxId = std::max(maxId, mId);
      }
      polyData->GetCellData()->AddArray(matIdArray);
      polyData->GetCellData()->SetActiveScalars("MaterialIds");
      VIENNACORE_LOG_DEBUG("Added MaterialIds array to cell data.");
    }

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    if (useMaterialIds) {
      mapper->SetScalarModeToUseCellData();
      mapper->ScalarVisibilityOn();
      mapper->SelectColorArray("MaterialIds");

      vtkSmartPointer<vtkLookupTable> lut =
          vtkSmartPointer<vtkLookupTable>::New();

      lut->SetNumberOfTableValues(256);
      lut->SetHueRange(0.667, 0.0); // blue → red
      lut->SetSaturationRange(1.0, 1.0);
      lut->SetValueRange(1.0, 1.0);
      lut->Build();

      mapper->SetLookupTable(lut);
      mapper->SetScalarRange(minId, maxId);
    }

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(3.0); // Thicker lines

    renderer->AddActor(actor);
    renderer->ResetCamera();
  }

  void updateVolumeMesh() {
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(volumeMesh);

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    renderer->AddActor(actor);
    renderer->ResetCamera();
  }

private:
  SmartPointer<Mesh<T>> mesh = nullptr;
  std::vector<T> materialIds;

  vtkSmartPointer<vtkUnstructuredGrid> volumeMesh = nullptr;

  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkRenderWindowInteractor> interactor;
  vtkSmartPointer<vtkPolyData> polyData;

  std::array<double, 3> backgroundColor = {84.0 / 255, 89.0 / 255, 109.0 / 255};
  std::array<int, 2> windowSize = {800, 600};
};

} // namespace viennals

#endif // VIENNALS_VTK_RENDERING