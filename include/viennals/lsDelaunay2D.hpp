#pragma once

#include <lsConstraintCleaner.hpp>
#include <lsDomain.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/write_VTU.h>

#include <vtkCellLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

namespace viennals {

using namespace viennacore;

template <typename NumericType> class Delaunay2D {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_vertex_base_2<K> Vb;
  typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

  SmartPointer<Mesh<NumericType>> mesh;
  std::vector<SmartPointer<Domain<NumericType, 2>>> domains;
  SmartPointer<MaterialMap> materialMap = nullptr;
  double maxTriangleSize = -1.;
  int bottomExtent = 1;
  int bottomLayerMaterialId = -1;
  bool closeDomain = true;
  bool cleanConstraints = true;
  bool verboseConstraintCleaning = false;
  NumericType constraintTargetSpacing = -1;
  NumericType constraintMergeThreshold = -1;
  NumericType constraintMinEdgeLength = -1;

private:
  void cdtToMesh(const CDT &cdt, const bool inDomain = true) {
    mesh->clear();
    // Use the address of the underlying CGAL vertex as a stable key.
    // (CGAL handles are not guaranteed to be hashable.)
    std::unordered_map<const void *, unsigned> vertexToPointId;
    vertexToPointId.reserve(static_cast<std::size_t>(cdt.number_of_vertices()));

    for (auto vit = cdt.finite_vertices_begin();
         vit != cdt.finite_vertices_end(); ++vit) {
      const auto &p = vit->point();
      const unsigned pid = mesh->insertNextNode(
          {CGAL::to_double(p.x()), CGAL::to_double(p.y()), 0.0});
      vertexToPointId.emplace(static_cast<const void *>(&*vit), pid);
    }

    CGAL::internal::In_domain<CDT> in_domain;

    for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(),
                                             end = cdt.finite_faces_end();
         fit != end; ++fit) {
      std::array<unsigned, 3> ids;
      if (!inDomain || get(in_domain, fit)) {
        bool ok = true;
        for (int i = 0; i < 3; ++i) {
          const auto vh = fit->vertex(i);
          const auto it = vertexToPointId.find(static_cast<const void *>(&*vh));
          if (it == vertexToPointId.end()) {
            ok = false;
            break;
          }
          ids[i] = it->second;
        }

        if (!ok)
          continue;

        mesh->insertNextTriangle(ids);
      }
    }
  }

  struct ExtremeIndices {
    unsigned minX_maxY;
    unsigned maxX_maxY;
  };

  ExtremeIndices
  findExtremePointIndices(const std::vector<Vec3D<NumericType>> &pts) {
    ExtremeIndices res{0, 0};
    if (pts.empty())
      return res;

    for (std::size_t i = 1; i < pts.size(); ++i) {
      const auto &p = pts[i];

      const auto &minP = pts[res.minX_maxY];
      if (p[0] < minP[0] || (p[0] == minP[0] && p[1] > minP[1])) {
        res.minX_maxY = i;
      }

      const auto &maxP = pts[res.maxX_maxY];
      if (p[0] > maxP[0] || (p[0] == maxP[0] && p[1] > maxP[1])) {
        res.maxX_maxY = i;
      }
    }

    return res;
  }

  void closeMesh() {
    // close mesh
    auto const gridDelta = domains.back()->getGrid().getGridDelta();
    auto const extremeIndices = findExtremePointIndices(mesh->nodes);
    auto const &minExtent = mesh->minimumExtent;
    auto const &maxExtent = mesh->maximumExtent;

    auto p1 = mesh->insertNextNode(
        {minExtent[0], minExtent[1] - bottomExtent * gridDelta, 0.});
    auto p2 = mesh->insertNextNode(
        {maxExtent[0], minExtent[1] - bottomExtent * gridDelta, 0.});

    mesh->insertNextLine({extremeIndices.maxX_maxY, p2});
    mesh->insertNextLine({p2, p1});
    mesh->insertNextLine({p1, extremeIndices.minX_maxY});
  }

  void createConstraints(CDT &cdt) {
    std::vector<CDT::Vertex_handle> vertexMap(mesh->nodes.size());
    for (unsigned i = 0; i < mesh->nodes.size(); ++i) {
      const auto &node = mesh->nodes[i];
      vertexMap[i] = cdt.insert(CDT::Point(node[0], node[1]));
    }

    auto const &lines = mesh->lines;
    for (const auto &line : lines) {
      cdt.insert_constraint(vertexMap[line[0]], vertexMap[line[1]]);
    }
  }

public:
  Delaunay2D() = default;

  void insertNextLevelSet(SmartPointer<Domain<NumericType, 2>> domain) {
    domains.push_back(domain);
  }

  void setMesh(SmartPointer<Mesh<NumericType>> passedMesh) {
    mesh = passedMesh;
  }

  void setMaxTriangleSize(double size) { maxTriangleSize = size; }

  void setBottomExtent(int extent) { bottomExtent = extent; }

  void setBottomLayerMaterialId(int materialId) {
    bottomLayerMaterialId = materialId;
  }

  void setMaterialMap(SmartPointer<MaterialMap> matMap) {
    materialMap = matMap;
  }

  void setCloseDomain(bool close) { closeDomain = close; }

  /// Enable/disable constraint cleaning before CDT
  void setCleanConstraints(bool clean) { cleanConstraints = clean; }

  /// Enable verbose output for constraint cleaning
  void setVerboseConstraintCleaning(bool verbose) {
    verboseConstraintCleaning = verbose;
  }

  /// Set target edge spacing for constraint cleaning (auto if < 0)
  void setConstraintTargetSpacing(NumericType spacing) {
    constraintTargetSpacing = spacing;
  }

  /// Set merge threshold for near-duplicate vertices (auto if < 0)
  void setConstraintMergeThreshold(NumericType threshold) {
    constraintMergeThreshold = threshold;
  }

  /// Set minimum edge length for constraint cleaning (auto if < 0)
  void setConstraintMinEdgeLength(NumericType length) {
    constraintMinEdgeLength = length;
  }

  void apply() {
    mesh->clear();

    ToMultiSurfaceMesh<NumericType, 2> mesher;
    WriteVisualizationMesh<NumericType, 2> visMesh;
    mesher.setMesh(mesh);
    for (const auto &d : domains) {
      mesher.insertNextLevelSet(d);
      visMesh.insertNextLevelSet(d);
    }
    visMesh.setWriteToFile(false);
    mesher.apply();
    visMesh.apply();

    // surface line mesh is now in mesh

    // remove normals
    mesh->getCellData().clear();

    // Clean constraints before CDT if enabled
    if (cleanConstraints) {
      ConstraintCleaner<NumericType> cleaner;
      cleaner.setPoints(mesh->nodes);
      cleaner.setEdges(mesh->lines);
      cleaner.setVerbose(verboseConstraintCleaning);

      if (constraintTargetSpacing > 0) {
        cleaner.setTargetSpacing(constraintTargetSpacing);
      }
      if (constraintMergeThreshold > 0) {
        cleaner.setMergeThreshold(constraintMergeThreshold);
      }
      if (constraintMinEdgeLength > 0) {
        cleaner.setMinEdgeLength(constraintMinEdgeLength);
      }

      cleaner.apply();
      cleaner.applyToMesh(mesh);
    }

    if (closeDomain)
      closeMesh();

    VTKWriter<NumericType>(mesh, "surface_mesh").apply();

    // create constraints
    CDT cdt;
    createConstraints(cdt);

    // run meshing
    auto const gridDelta = domains.back()->getGrid().getGridDelta();
    maxTriangleSize = std::max(maxTriangleSize, gridDelta);
    CGAL::refine_Delaunay_mesh_2(
        cdt, CGAL::parameters::criteria(Criteria(0.125, maxTriangleSize)));

    // std::fstream ofs("cdt_mesh.vtu", std::ios::out);
    // CGAL::IO::write_VTU(ofs, cdt);

    auto rgrid = visMesh.getVolumeMesh();
    auto materials = rgrid->GetCellData()->GetArray("Material");

    auto cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(rgrid);
    cellLocator->BuildLocator();

    cdtToMesh(cdt);

    std::vector<double> materialIds;
    for (const auto &tri : mesh->triangles) {
      auto &p1 = mesh->nodes[tri[0]];
      auto &p2 = mesh->nodes[tri[1]];
      auto &p3 = mesh->nodes[tri[2]];
      std::array<double, 3> centroid = {(p1[0] + p2[0] + p3[0]) / 3.,
                                        (p1[1] + p2[1] + p3[1]) / 3.,
                                        (p1[2] + p2[2] + p3[2]) / 3.};

      vtkIdType cellId = cellLocator->FindCell(centroid.data());

      if (cellId == -1) {
        if (bottomLayerMaterialId == -1) {
          double materialId = 0.0;
          if (materialMap) {
            materialId = materialMap->getMaterialId(0);
          }
          materialIds.push_back(materialId);
        } else {
          materialIds.push_back(static_cast<double>(bottomLayerMaterialId));
        }
      } else {
        double materialId = materials->GetTuple1(cellId);
        if (materialMap) {
          materialId =
              materialMap->getMaterialId(static_cast<std::size_t>(materialId));
        }
        materialIds.push_back(materialId);
      }
    }
    mesh->getCellData().insertNextScalarData(materialIds, "Material");
  }
};
} // namespace viennals