#include <CreateSurfaceMesh.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/write_VTU.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const double radius = 15.0;
  const ls::VectorType<double, D> centre{0., 0.};

  ls::MakeGeometry<double, 2>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
      .apply();

  viennaps::CreateSurfaceMesh<double, double, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "sphere.vtp").apply();

  CDT cdt;

  std::vector<Vertex_handle> vertexMap(mesh->nodes.size());
  for (unsigned i = 0; i < mesh->nodes.size(); ++i) {
    const auto &node = mesh->nodes[i];
    vertexMap[i] = cdt.insert(Point(node[0], node[1]));
  }

  auto const &lines = mesh->lines;
  for (const auto &line : lines) {
    cdt.insert_constraint(vertexMap[line[0]], vertexMap[line[1]]);
  }

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the triangulation..." << std::endl;
  CGAL::refine_Delaunay_mesh_2(cdt,
                               CGAL::parameters::criteria(Criteria(0.125, 2.)));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::fstream ofs("cdt_mesh.vtu", std::ios::out);
  CGAL::IO::write_VTU(ofs, cdt);

  return 0;
}
