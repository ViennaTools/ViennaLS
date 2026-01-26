#pragma once

#include <lsDomain.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

namespace viennals {

using namespace viennacore;

template <typename NumericType> class Delaunay3D {

  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;

  using Vbb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
  using Vb =
      CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vbb>;

  using Cbb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
  using Cb =
      CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cbb>;

  using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Tr = CGAL::Triangulation_3<K, Tds>;
  using CCDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;

  // Triangulation for Remeshing
  using CCDT_Tr = CCDT::Triangulation;
  using Triangulation_3 =
      CGAL::Triangulation_3<K, CCDT_Tr::Triangulation_data_structure>;

  using Vertex_handle = Triangulation_3::Vertex_handle;
  using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
  using Constraints_set =
      std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
  using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

  SmartPointer<Mesh<NumericType>> mesh;
  std::vector<SmartPointer<Domain<NumericType, 3>>> domains;
  SmartPointer<MaterialMap> materialMap = nullptr;

public:
  Delaunay3D() = default;

  Delaunay3D(SmartPointer<Mesh<NumericType>> passedMesh) : mesh(passedMesh) {}

  void insertNextLevelSet(SmartPointer<Domain<NumericType, 3>> domain) {
    domains.push_back(domain);
  }

  void setMesh(SmartPointer<Mesh<NumericType>> passedMesh) {
    mesh = passedMesh;
  }

  void apply() {

    ToMultiSurfaceMesh<NumericType, 3> converter;
    WriteVisualizationMesh<NumericType, 3> visMesh;
    for (const auto &domain : domains) {
      visMesh.insertNextLevelSet(domain);
      converter.insertNextLevelSet(domain);
    }
    visMesh.setExtractHullMesh(true);
    visMesh.setExtractVolumeMesh(true);
    visMesh.setFileName("delaunay3D_initial_mesh");
    visMesh.apply();

    converter.setMesh(mesh);
    converter.apply();

    VTKWriter<NumericType>(mesh, "delaunay3D_surface_mesh").apply();

    auto hullVTK = visMesh.getHullMesh();
    auto materialIds = hullVTK->GetCellData()->GetArray("Material");

    {
      std::vector<std::array<unsigned, 3>> triangles;
      std::vector<Vec3D<NumericType>> points;
      for (vtkIdType i = 0; i < hullVTK->GetNumberOfPoints(); ++i) {
        double p[3];
        hullVTK->GetPoint(i, p);
        points.push_back(Vec3D<NumericType>{static_cast<NumericType>(p[0]),
                                            static_cast<NumericType>(p[1]),
                                            static_cast<NumericType>(p[2])});
      }

      for (vtkIdType i = 0; i < hullVTK->GetNumberOfCells(); ++i) {
        vtkCell *cell = hullVTK->GetCell(i);
        int materialId = static_cast<int>(materialIds->GetTuple1(i));
        if (cell->GetNumberOfPoints() == 3 && materialId == 1) {
          std::array<unsigned, 3> tri;
          for (vtkIdType j = 0; j < 3; ++j) {
            tri[j] = static_cast<unsigned>(cell->GetPointId(j));
          }
          triangles.push_back(tri);
        }
      }
      mesh->nodes = std::move(points);
      mesh->triangles = std::move(triangles);

      VTKWriter<NumericType>(mesh, "delaunay3D_surface_mesh").apply();
    }

    // std::vector<K::Point_3> points;
    // std::vector<std::vector<std::size_t>> polygons;
    // points.reserve(mesh->getNodes().size());
    // polygons.reserve(mesh->triangles.size());
    // for (const auto &node : mesh->getNodes()) {
    //   points.emplace_back(node[0], node[1], node[2]);
    // }
    // for (const auto &tri : mesh->triangles) {
    //   polygons.push_back({tri[0], tri[1], tri[2]});
    // }

    SurfaceMesh cgalMesh;
    std::vector<typename SurfaceMesh::Vertex_index> vertexHandles;
    vertexHandles.reserve(mesh->nodes.size());
    for (const auto &point : mesh->nodes) {
      vertexHandles.push_back(
          cgalMesh.add_vertex(K::Point_3(point[0], point[1], point[2])));
    }
    for (const auto &tri : mesh->triangles) {
      cgalMesh.add_face(vertexHandles[tri[0]], vertexHandles[tri[1]],
                        vertexHandles[tri[2]]);
    }

    // unsigned maxX_maxY_maxZ = 0;
    // unsigned maxX_maxY_minZ = 0;
    // unsigned maxX_minY_maxZ = 0;
    // unsigned maxX_minY_minZ = 0;
    // unsigned minX_maxY_maxZ = 0;
    // unsigned minX_maxY_minZ = 0;
    // unsigned minX_minY_maxZ = 0;
    // unsigned minX_minY_minZ = 0;

    // for (auto const &node : mesh->getNodes()) {
    //   points.emplace_back(node[0], node[1], node[2]);

    //   if (node[0] >= points[maxX_maxY_maxZ].x() &&
    //       node[1] >= points[maxX_maxY_maxZ].y() &&
    //       node[2] >= points[maxX_maxY_maxZ].z()) {
    //     maxX_maxY_maxZ = points.size() - 1;
    //   }
    //   if (node[0] >= points[maxX_maxY_minZ].x() &&
    //       node[1] >= points[maxX_maxY_minZ].y() &&
    //       node[2] <= points[maxX_maxY_minZ].z()) {
    //     maxX_maxY_minZ = points.size() - 1;
    //   }
    //   if (node[0] >= points[maxX_minY_maxZ].x() &&
    //       node[1] <= points[maxX_minY_maxZ].y() &&
    //       node[2] >= points[maxX_minY_maxZ].z()) {
    //     maxX_minY_maxZ = points.size() - 1;
    //   }
    //   if (node[0] >= points[maxX_minY_minZ].x() &&
    //       node[1] <= points[maxX_minY_minZ].y() &&
    //       node[2] <= points[maxX_minY_minZ].z()) {
    //     maxX_minY_minZ = points.size() - 1;
    //   }
    //   if (node[0] <= points[minX_maxY_maxZ].x() &&
    //       node[1] >= points[minX_maxY_maxZ].y() &&
    //       node[2] >= points[minX_maxY_maxZ].z()) {
    //     minX_maxY_maxZ = points.size() - 1;
    //   }
    //   if (node[0] <= points[minX_maxY_minZ].x() &&
    //       node[1] >= points[minX_maxY_minZ].y() &&
    //       node[2] <= points[minX_maxY_minZ].z()) {
    //     minX_maxY_minZ = points.size() - 1;
    //   }
    //   if (node[0] <= points[minX_minY_maxZ].x() &&
    //       node[1] <= points[minX_minY_maxZ].y() &&
    //       node[2] >= points[minX_minY_maxZ].z()) {
    //     minX_minY_maxZ = points.size() - 1;
    //   }
    //   if (node[0] <= points[minX_minY_minZ].x() &&
    //       node[1] <= points[minX_minY_minZ].y() &&
    //       node[2] <= points[minX_minY_minZ].z()) {
    //     minX_minY_minZ = points.size() - 1;
    //   }

    //   for (const auto &tri : mesh->triangles) {
    //     polygons.push_back({tri[0], tri[1], tri[2]});
    //   }
    // }

    // namespace PMP = CGAL::Polygon_mesh_processing;
    // PMP::polygon_mesh_to_polygon_soup(cgalMesh, points, polygons);
    // PMP::autorefine_triangle_soup(points, polygons);

    std::cout << "Starting Delaunay tetrahedralization..." << std::endl;
    CCDT ccdt =
        CGAL::make_conforming_constrained_Delaunay_triangulation_3<CCDT>(
            cgalMesh);

    double target_edge_length = domains[0]->getGrid().getGridDelta() * 2.0;
    unsigned int iterations = 5;

    std::cout << "Number of vertices in the CDT: "
              << ccdt.triangulation().number_of_vertices() << '\n';
    std::cout << "Number of constrained facets in the CDT: "
              << ccdt.number_of_constrained_facets() << '\n';

    Constraints_set constraints;
    Constraints_pmap constraints_pmap(constraints);

    namespace np = CGAL::parameters;
    namespace Tet_remesh = CGAL::Tetrahedral_remeshing;
    Tr tr = Tet_remesh::get_remeshing_triangulation(
        std::move(ccdt), np::edge_is_constrained_map(constraints_pmap));

    std::cout << "Starting isotropic remeshing..." << std::endl;
    CGAL::tetrahedral_isotropic_remeshing(
        tr, target_edge_length,
        np::number_of_iterations(iterations)
            .edge_is_constrained_map(constraints_pmap));

    std::cout << "There are " << tr.number_of_vertices()
              << " vertices after remeshing" << std::endl;

    mesh->clear();
    for (auto vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit) {
      const auto &p = vit->point();
      mesh->insertNextNode({CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                            CGAL::to_double(p.z())});
    }

    std::vector<std::array<unsigned, 4>> tetras;

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end();
         ++cit) {
      // if cell not in domain, skip
      bool cell_in_domain = false;
      for (int i = 0; i < 4; ++i) {
        auto fh = std::make_pair(cit, i);
        if (constraints.find(fh) != constraints.end()) {
          cell_in_domain = true;
          break;
        }
      }
      if (!cell_in_domain)
        continue;

      std::array<unsigned, 4> tet;
      for (int i = 0; i < 4; ++i) {
        tet[i] = static_cast<unsigned>(
            std::distance(tr.vertices_begin(), cit->vertex(i)));
      }
      tetras.push_back(tet);
    }

    // using Tr_ = std::remove_cv_t<std::remove_reference_t<decltype(tr)>>;
    // using Vertex_handle = typename Tr_::Vertex_handle;
    // using Facet = typename Tr_::Facet;
    // using Cell_handle = typename Tr_::Cell_handle;
    // CGAL::unordered_flat_map<Cell_handle, bool> cells_in_domain;
    // for (Cell_handle c : tr.all_cell_handles())
    //   cells_in_domain[c] = true;

    // std::stack<Cell_handle> stack;
    // stack.push(tr.infinite_cell());
    // while (!stack.empty()) {
    //   auto ch = stack.top();
    //   stack.pop();
    //   cells_in_domain[ch] = false;
    //   for (int i = 0; i < 4; ++i) {
    //     if (ccdt.is_facet_constrained(ch, i))
    //       continue;
    //     auto n = ch->neighbor(i);
    //     if (cells_in_domain[n])
    //       stack.push(n);
    //   }
    // }

    // for (auto &cell : cells_in_domain) {
    //   if (cell.second) {
    //     std::array<unsigned, 4> tet;
    //     for (int i = 0; i < 4; ++i) {
    //       tet[i] = static_cast<unsigned>(
    //           std::distance(tr.vertices_begin(), cell.first->vertex(i)));
    //     }
    //     tetras.push_back(tet);
    //   }
    // }

    std::cout << "Number of tets in the domain: " << tetras.size() << std::endl;
    mesh->tetras = std::move(tetras);
    VTKWriter<NumericType>(mesh, "delaunay3D_mesh.vtu").apply();

    auto rgrid = visMesh.getVolumeMesh();
    auto materials = rgrid->GetCellData()->GetArray("Material");
  }
};
} // namespace viennals