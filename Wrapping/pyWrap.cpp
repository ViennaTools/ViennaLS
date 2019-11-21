/*
  This file is used to generate the python module of viennals.
  It uses pybind11 to create the modules.

  All necessary headers are included here and the interface
  of the classes which should be exposed defined
*/


#include <pybind11/pybind11.h>

// all header files which define API functions
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCheck.hpp>
#include <lsConvexHull.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

typedef double T;
constexpr int D = 2;

// module specification
PYBIND11_MODULE(viennaLS, module) {
  module.doc() = "ViennaLS python module.";

  // // lsAdvect
  // pybind11::class_< lsAdvect<double, 2> >(module, "lsAdvect")
  //     // constructors
  //     .def(pybind11::init<>())
  //     // getters and setters
  //     .def("getAdvectionTime", &lsAdvect<double, 2>::getAdvectionTime, "Get the advection time.");
  //     .def("setAdvectionTime", &lsAdvect<double, 2>::setAdvectionTime, "Set the time until when the level set should be advected.");
  //     // member methods
  //     // .def("getName", (std::string &(Pet::*)()) &Pet::getName, "Get a mutable reference to name.")
  //     // .def("getName", (const std::string &(Pet::*)() const) &Pet::getName, "Get const reference.");

  // lsDomain
  pybind11::class_< lsDomain<T, D> >(module, "lsDomain")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<hrleCoordType>())
      .def(pybind11::init<hrleCoordType *, lsDomain<T, D>::BoundaryType *>())
      .def(pybind11::init<hrleCoordType *, lsDomain<T, D>::BoundaryType *, hrleCoordType>())
      .def(pybind11::init<lsDomain<T, D>::PointValueVectorType,  hrleCoordType *,
               lsDomain<T, D>::BoundaryType *>())
      .def(pybind11::init<lsDomain<T, D>::PointValueVectorType,  hrleCoordType *,
               lsDomain<T, D>::BoundaryType *, hrleCoordType>())

      .def(pybind11::init<const lsDomain<T, D> &>())

      // methods
      .def("print", &lsDomain<T, D>::print, "Print level set structure.")
      ;

  // lsSphere
  pybind11::class_< lsSphere<T, D> >(module, "lsSphere")
      // constructors
      .def(pybind11::init<T, T, T>())
      .def(pybind11::init<T, T, T, T>())
      // methods
      ;

  // lsMesh
  pybind11::class_< lsMesh >(module, "lsMesh")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("print", &lsMesh::print, "Print basic information about the mesh.")
      ;

  // lsMakeGeometry
  pybind11::class_< lsMakeGeometry<T, D> >(module, "lsMakeGeometry")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsSphere<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsPlane<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsBox<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsPointCloud<T, D> &>())
      // methods
      .def("setLevelSet", &lsMakeGeometry<T, D>::setLevelSet, "Set the levelset in which to create the geometry.")
      .def("setGeometry", (void (lsMakeGeometry<T, D>::*)(lsSphere<T, D>&)) &lsMakeGeometry<T, D>::setGeometry)

      .def("apply", &lsMakeGeometry<T, D>::apply, "Generate the geometry.")
      ;

  // lsToSurfaceMesh
  pybind11::class_< lsToSurfaceMesh<T, D> >(module, "lsToSurfaceMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&>())
      // methods
      .def("apply", &lsToSurfaceMesh<T, D>::apply, "Convert the levelset to a surface mesh.")
      ;

  // lsVTKWriter
  pybind11::class_< lsVTKWriter >(module, "lsVTKWriter")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh&>())
      .def(pybind11::init<lsMesh&, std::string>())
      // methods
      .def("apply", &lsVTKWriter::apply, "Write the mesh.")
      ;

}
