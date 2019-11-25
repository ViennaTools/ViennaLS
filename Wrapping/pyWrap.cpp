/*
  This file is used to generate the python module of viennals.
  It uses pybind11 to create the modules.

  All necessary headers are included here and the interface
  of the classes which should be exposed defined
*/

// correct module name macro
#define TOKENPASTE_INTERNAL(x, y, z) x ## y ## z
#define TOKENPASTE(x, y, z) TOKENPASTE_INTERNAL(x, y, z)
#define VIENNALS_MODULE_NAME TOKENPASTE(viennaLS, VIENNALS_PYTHON_DIMENSION, d)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// all header files which define API functions
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCheck.hpp>
#include <lsConvexHull.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFileFormats.hpp>
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

// always use double for python export
typedef double T;
// get dimension from cmake define
constexpr int D = VIENNALS_PYTHON_DIMENSION;

// define trampoline classes for interface functions
// ALSO NEED TO ADD TRAMPOLINE CLASSES FOR CLASSES
// WHICH HOLD REFERENCES TO INTERFACE(ABSTRACT) CLASSES

// BASE CLASS WRAPPERS
// lsVelocityField only defines interface and has no functionality
class PylsVelocityField : public lsVelocityField<T> {
  typedef hrleVectorType<T, 3> vectorType;
  using lsVelocityField<T>::lsVelocityField;
public:
  T getScalarVelocity(
      vectorType coordinate, int material,
      vectorType normalVector = vectorType(T(0))) override {
        PYBIND11_OVERLOAD_PURE(T, lsVelocityField<T>, getScalarVelocity, coordinate, material, normalVector);
      }

  vectorType getVectorVelocity(
      hrleVectorType<T, 3> coordinate, int material,
      hrleVectorType<T, 3> normalVector = hrleVectorType<T, 3>(T(0))) override {
        PYBIND11_OVERLOAD_PURE(vectorType, lsVelocityField<T>, getVectorVelocity, coordinate, material, normalVector);
      }
};

//REFERNCE HOLDING CLASS WRAPPERS

// maybe needed wrapper code once we move to smart pointers
// https://github.com/pybind/pybind11/issues/1389
// // lsAdvect wrapping since it holds lsVelocityField references
// class PylsAdvect : public lsAdvect<T, D> {
//   pybind11::object pyObj;
// public:
//   PylsAdvect(lsDomain<T, D> &passedDomain, lsVelocityField<T> &passedVelocities) : lsAdvect<T,D>(passedDomain, passedVelocities), pyObj(pybind11::cast(passedVelocities)) {}
//
//   PylsAdvect(lsVelocityField<T> &passedVelocities) : lsAdvect<T,D>(passedVelocities), pyObj(pybind11::cast(passedVelocities)) {}
// };

// module specification
PYBIND11_MODULE(VIENNALS_MODULE_NAME, module) {
  module.doc() = "ViennaLS python module.";

  // maximum number of threads to use
  // must be set to one for now because lsAdvect
  // hangs a the end of parallel section
  // most likely to do with some omp issue
  omp_set_num_threads(1);

  // hrleVectorType
  pybind11::class_< hrleVectorType<T, 3> >(module, "lsVectorType")
      // constructors
      .def(pybind11::init<double, double, double>())
      ;

  // lsAdvect
  pybind11::class_< lsAdvect<T, D> >(module, "lsAdvect")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsVelocityField<T> &>())
      .def(pybind11::init<lsVelocityField<T> &>())
      // getters and setters
      .def("insertNextLevelSet", &lsAdvect<T, D>::insertNextLevelSet, "Insert next level set to use for advection.")
      .def("setVelocityField", &lsAdvect<T, D>::setVelocityField, "Set the velocity to use for advection.")
      .def("setAdvectionTime", &lsAdvect<T, D>::setAdvectionTime, "Set the time until when the level set should be advected.")
      .def("setTimeStepRatio", &lsAdvect<T, D>::setTimeStepRatio, "Set the maximum time step size relative to grid size. Advection is only stable for <0.5.")
      .def("setCalculateNormalVectors", &lsAdvect<T, D>::setCalculateNormalVectors, "Set whether normal vectors are needed for the supplied velocity field.")
      .def("setIgnoreVoids", &lsAdvect<T, D>::setIgnoreVoids, "Set whether voids in the geometry should be ignored during advection or not.")
      .def("getAdvectionTime", &lsAdvect<T, D>::getAdvectionTime, "Get the advection time.")
      .def("getNumberOfTimeSteps", &lsAdvect<T, D>::getNumberOfTimeSteps, "Get how many advection steps were performed after the last apply() call.")
      .def("getTimeStepRatio", &lsAdvect<T, D>::getTimeStepRatio, "Get the time step ratio used for advection.")
      .def("getCalculateNormalVectors", &lsAdvect<T, D>::getCalculateNormalVectors, "Get whether normal vectors are computed during advection.")
      .def("setIntegrationScheme", &lsAdvect<T, D>::setIntegrationScheme, "Set the integration scheme to use during advection.")
      .def("setDissipationAlpha", &lsAdvect<T, D>::setDissipationAlpha, "Set the dissipation value to use for Lax Friedrichs integration.")
      .def("apply", &lsAdvect<T, D>::apply, "Perform advection.")
      ;
  // enums
  pybind11::enum_<lsIntegrationSchemeEnum>(module, "lsIntegrationSchemeEnum")
      .value("ENGQUIST_OSHER_1ST_ORDER", lsIntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER)
      .value("ENGQUIST_OSHER_2ND_ORDER", lsIntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER)
      .value("LAX_FRIEDRICHS_1ST_ORDER", lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER)
      .value("LAX_FRIEDRICHS_2ND_ORDER", lsIntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER)
      .value("STENCIL_LOCAL_LAX_FRIEDRICHS", lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS)
      ;

  // lsBooleanOperation
  pybind11::class_< lsBooleanOperation<T, D> >(module, "lsBooleanOperation")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsBooleanOperationEnum>())
      .def(pybind11::init<lsDomain<T, D> &, lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsDomain<T, D> &, lsBooleanOperationEnum>())
      // methods
      .def("setLevelset", &lsBooleanOperation<T, D>::setLevelSet, "Set levelset on which the boolean operation should be performed.")
      .def("setSecondLevelSet", &lsBooleanOperation<T, D>::setSecondLevelSet, "Set second levelset for boolean operation.")
      .def("setBooleanOperation", &lsBooleanOperation<T, D>::setBooleanOperation, "Set which type of boolean operation should be performed.")
      .def("apply", &lsBooleanOperation<T, D>::apply, "Perform the boolean operation.")
      ;
  // enums
  pybind11::enum_<lsBooleanOperationEnum>(module, "lsBooleanOperationEnum")
      .value("INTERSECT", lsBooleanOperationEnum::INTERSECT)
      .value("UNION", lsBooleanOperationEnum::UNION)
      .value("RELATIVE_COMPLEMENT", lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .value("INVERT", lsBooleanOperationEnum::INVERT)
      ;

  // lsCalculateNormalVectors
  pybind11::class_< lsCalculateNormalVectors<T, D> >(module, "lsCalculateNormalVectors")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, bool>())
      //methods
      .def("setLevelSet", &lsCalculateNormalVectors<T, D>::setLevelSet, "Set levelset for which to calculate normal vectors.")
      .def("setOnlyActivePoints", &lsCalculateNormalVectors<T, D>::setOnlyActivePoints, "Set whether normal vectors should only be calculated for level set points <0.5.")
      .def("apply", &lsCalculateNormalVectors<T, D>::apply, "Perform normal vector calculation.")
      ;

  // lsCheck
  pybind11::class_< lsCheck<T, D> >(module, "lsCheck")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      //methods
      .def("setLevelSet", &lsCheck<T, D>::setLevelSet, "Set levelset for which to calculate normal vectors.")
      .def("apply", &lsCheck<T, D>::apply, "Perform check.")
      ;

  // lsConvexHull
  pybind11::class_< lsConvexHull<T, D> >(module, "lsConvexHull")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh &, lsPointCloud<T, D>&>())
      //methods
      .def("setMesh", &lsConvexHull<T, D>::setMesh, "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &lsConvexHull<T, D>::setPointCloud, "Set point cloud used to generate mesh.")
      .def("apply", &lsConvexHull<T, D>::apply, "Perform check.")
      ;

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

  // lsExpand
  pybind11::class_< lsExpand<T, D> >(module, "lsExpand")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, int>())
      //methods
      .def("setLevelSet", &lsExpand<T, D>::setLevelSet, "Set levelset to expand.")
      .def("setWidth", &lsExpand<T, D>::setWidth, "Set the width to expand to.")
      .def("apply", &lsExpand<T, D>::apply, "Perform expansion.")
      ;

  // lsFileFormats
  pybind11::enum_<lsFileFormatEnum>(module, "lsFileFormatEnum")
      .value("VTK_LEGACY", lsFileFormatEnum::VTK_LEGACY)
      .value("VTP", lsFileFormatEnum::VTP)
      .value("VTU", lsFileFormatEnum::VTU)
      ;

  // lsFromSurfaceMesh
  pybind11::class_< lsFromSurfaceMesh<T, D> >(module, "lsFromSurfaceMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D>&, lsMesh&>())
      .def(pybind11::init<lsDomain<T, D>&, lsMesh&, bool>())
      // methods
      .def("setLevelSet", &lsFromSurfaceMesh<T, D>::setLevelSet, "Set levelset to read into.")
      .def("setMesh", &lsFromSurfaceMesh<T, D>::setMesh, "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles", &lsFromSurfaceMesh<T, D>::setRemoveBoundaryTriangles, "Set whether to include mesh elements outside of the simulation domain.")
      .def("apply", &lsFromSurfaceMesh<T, D>::apply, "Construct a levelset from a surface mesh.")
      ;

  // lsFromVolumeMesh
  pybind11::class_< lsFromVolumeMesh<T, D> >(module, "lsFromVolumeMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<std::vector<lsDomain<T, D>> &, lsMesh &>())
      .def(pybind11::init<std::vector<lsDomain<T, D>> &, lsMesh &, bool>())
      // methods
      .def("setLevelSets", &lsFromVolumeMesh<T, D>::setLevelSets, "Set levelsets to read into.")
      .def("setMesh", &lsFromVolumeMesh<T, D>::setMesh, "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles", &lsFromVolumeMesh<T, D>::setRemoveBoundaryTriangles, "Set whether to include mesh elements outside of the simulation domain.")
      .def("apply", &lsFromVolumeMesh<T, D>::apply, "Construct a levelset from a volume mesh.")
      ;

  // lsGeometries
  // lsSphere
  pybind11::class_< lsSphere<T, D> >(module, "lsSphere")
      // constructors
      .def(pybind11::init<const std::vector<T> &, T>())
      // methods
      ;
  // lsPlane
  pybind11::class_< lsPlane<T, D> >(module, "lsPlane")
      // constructors
      .def(pybind11::init<const std::vector<T> &, const std::vector<T> &>())
      // methods
      ;

  // lsBox
  pybind11::class_< lsBox<T, D> >(module, "lsBox")
      // constructors
      .def(pybind11::init<const std::vector<T> &, const std::vector<T> &>())
      // methods
      ;

  // lsPointCloud
  pybind11::class_< lsPointCloud<T, D> >(module, "lsPointCloud")
      // constructors
      .def(pybind11::init<const std::vector<std::vector<T> > &>())
      // methods
      .def("insertNextPoint", (void (lsPointCloud<T,D>::*)(const std::vector<T> &)) &lsPointCloud<T, D>::insertNextPoint);

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

  // lsMesh
  pybind11::class_< lsMesh >(module, "lsMesh")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("print", &lsMesh::print, "Print basic information about the mesh.")
      ;

  // lsPrune
  pybind11::class_< lsPrune<T, D> >(module, "lsPrune")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      //methods
      .def("setLevelSet", &lsPrune<T, D>::setLevelSet, "Set levelset to prune.")
      .def("apply", &lsPrune<T, D>::apply, "Perform pruning operation.")
      ;

  // lsReduce
  pybind11::class_< lsReduce<T, D> >(module, "lsReduce")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, int>())
      .def(pybind11::init<lsDomain<T, D> &, int, bool>())
      //methods
      .def("setLevelSet", &lsReduce<T, D>::setLevelSet, "Set levelset to reduce.")
      .def("setWidth", &lsReduce<T, D>::setWidth, "Set the width to reduce to.")
      .def("setNoNewSegment", &lsReduce<T, D>::setNoNewSegment, "Set whether the levelset should be segmented anew (balanced across cores) after reduction.")
      .def("apply", &lsReduce<T, D>::apply, "Perform reduction.")
      ;

  // lsToDiskMesh
  pybind11::class_< lsToDiskMesh<T, D> >(module, "lsToDiskMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D>&, lsMesh&>())
      // methods
      .def("setLevelSet", &lsToDiskMesh<T, D>::setLevelSet, "Set levelset to mesh.")
      .def("setMesh", &lsToDiskMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &lsToDiskMesh<T, D>::apply, "Convert the levelset to a surface mesh.")
      ;

  // lsToMesh
  pybind11::class_< lsToMesh<T, D> >(module, "lsToMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&, bool>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&, bool, bool>())
      // methods
      .def("setLevelSet", &lsToMesh<T, D>::setLevelSet, "Set levelset to mesh.")
      .def("setMesh", &lsToMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("setOnlyDefined", &lsToMesh<T, D>::setOnlyDefined, "Set whether only defined points should be output to the mesh.")
      .def("setOnlyActive", &lsToMesh<T, D>::setOnlyActive, "Set whether only level set points <0.5 should be output.")
      .def("apply", &lsToMesh<T, D>::apply, "Convert the levelset to a surface mesh.")
      ;

  // lsToSurfaceMesh
  pybind11::class_< lsToSurfaceMesh<T, D> >(module, "lsToSurfaceMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&>())
      // methods
      .def("setLevelSet", &lsToSurfaceMesh<T, D>::setLevelSet, "Set levelset to mesh.")
      .def("setMesh", &lsToSurfaceMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &lsToSurfaceMesh<T, D>::apply, "Convert the levelset to a surface mesh.")
      ;

  // lsToVoxelMesh
  pybind11::class_< lsToVoxelMesh<T, D> >(module, "lsToVoxelMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh&>())
      .def(pybind11::init<const lsDomain<T, D>&, lsMesh&>())
      .def(pybind11::init<const std::vector<const lsDomain<T, D> *> &, lsMesh&>())
      // methods
      .def("insertNextLevelSet", &lsToVoxelMesh<T, D>::insertNextLevelSet, "Insert next level set to output in the mesh.")
      .def("setMesh", &lsToVoxelMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &lsToVoxelMesh<T, D>::apply, "Convert the levelset to a surface mesh.")
      ;

  // lsVelocityField
  pybind11::class_< lsVelocityField<T>, PylsVelocityField >(module, "lsVelocityField")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("getScalarVelocity", &lsVelocityField<T>::getScalarVelocity, "Return the scalar velocity for a point of material at coordinate with normal vector normal.")
      .def("getVectorVelocity", &lsVelocityField<T>::getVectorVelocity, "Return the vector velocity for a point of material at coordinate with normal vector normal.")
      ;

  // lsVTKReader
  pybind11::class_< lsVTKReader >(module, "lsVTKReader")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh&>())
      .def(pybind11::init<lsMesh&, lsFileFormatEnum, std::string>())
      .def(pybind11::init<lsMesh&, std::string>())
      // methods
      .def("setMesh", &lsVTKReader::setMesh, "Set the mesh to read into.")
      .def("apply", &lsVTKReader::apply, "Read the mesh.")
      ;

  // lsVTKWriter
  pybind11::class_< lsVTKWriter >(module, "lsVTKWriter")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh&>())
      .def(pybind11::init<lsMesh&, lsFileFormatEnum, std::string>())
      .def(pybind11::init<lsMesh&, std::string>())
      // methods
      .def("setMesh", &lsVTKWriter::setMesh, "Set the mesh to output.")
      .def("apply", &lsVTKWriter::apply, "Write the mesh.")
      ;

}
