/*
  This file is used to generate the python module of viennals.
  It uses pybind11 to create the modules.

  All necessary headers are included here and the interface
  of the classes which should be exposed defined
*/

// correct module name macro
#define TOKENPASTE_INTERNAL(x, y, z) x##y##z
#define TOKENPASTE(x, y, z) TOKENPASTE_INTERNAL(x, y, z)
#define VIENNALS_MODULE_NAME TOKENPASTE(viennaLS, VIENNALS_PYTHON_DIMENSION, d)
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)
#define VIENNALS_MODULE_VERSION STRINGIZE(VIENNALS_VERSION)

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
#include <lsFastAdvect.hpp>
#include <lsFastAdvectDistributions.hpp>
#include <lsFileFormats.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsPointData.hpp>
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
  typedef std::array<T, 3> vectorType;
  using lsVelocityField<T>::lsVelocityField;

public:
  T getScalarVelocity(const vectorType &coordinate, int material,
                      const vectorType &normalVector) override {
    PYBIND11_OVERLOAD(T, lsVelocityField<T>, getScalarVelocity, coordinate,
                      material, normalVector);
  }

  vectorType getVectorVelocity(const vectorType &coordinate, int material,
                               const vectorType &normalVector) override {
    PYBIND11_OVERLOAD(vectorType, lsVelocityField<T>, getVectorVelocity,
                      coordinate, material, normalVector);
  }
};

// lsFastAdvectDistribution
class PylsFastAdvectDistribution : public lsFastAdvectDistribution<T, D> {
  typedef std::array<hrleCoordType, D> vectorType;
  typedef lsFastAdvectDistribution<T, D> ClassType;
  using lsFastAdvectDistribution<T, D>::lsFastAdvectDistribution;

public:
  bool isInside(const vectorType &v, double eps = 0.) const override {
    PYBIND11_OVERLOAD_PURE(bool, ClassType, isInside, v, eps);
  }

  T getSignedDistance(const vectorType &v) const override {
    PYBIND11_OVERLOAD_PURE(T, ClassType, getSignedDistance, v);
  }

  void getBounds(std::array<hrleCoordType, 2 * D> &bounds) const override {
    PYBIND11_OVERLOAD_PURE(void, ClassType, getBounds, bounds);
  }
};

// REFERENCE HOLDING CLASS WRAPPERS

// maybe needed wrapper code once we move to smart pointers
// https://github.com/pybind/pybind11/issues/1389
// // lsAdvect wrapping since it holds lsVelocityField references
// class PylsAdvect : public lsAdvect<T, D> {
//   pybind11::object pyObj;
// public:
//   PylsAdvect(lsDomain<T, D> &passedDomain, lsVelocityField<T>
//   &passedVelocities) : lsAdvect<T,D>(passedDomain, passedVelocities),
//   pyObj(pybind11::cast(passedVelocities)) {}
//
//   PylsAdvect(lsVelocityField<T> &passedVelocities) :
//   lsAdvect<T,D>(passedVelocities), pyObj(pybind11::cast(passedVelocities)) {}
// };

// module specification
PYBIND11_MODULE(VIENNALS_MODULE_NAME, module) {
  module.doc() = "ViennaLS python module.";

  // set version string of python module
  module.attr("__version__") = VIENNALS_MODULE_VERSION;

  // lsAdvect
  pybind11::class_<lsAdvect<T, D>>(module, "lsAdvect")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsVelocityField<T> &>())
      .def(pybind11::init<lsVelocityField<T> &>())
      // getters and setters
      .def("insertNextLevelSet", &lsAdvect<T, D>::insertNextLevelSet,
           "Insert next level set to use for advection.")
      .def("setVelocityField", &lsAdvect<T, D>::setVelocityField,
           "Set the velocity to use for advection.")
      .def("setAdvectionTime", &lsAdvect<T, D>::setAdvectionTime,
           "Set the time until when the level set should be advected.")
      .def("setTimeStepRatio", &lsAdvect<T, D>::setTimeStepRatio,
           "Set the maximum time step size relative to grid size. Advection is "
           "only stable for <0.5.")
      .def("setCalculateNormalVectors",
           &lsAdvect<T, D>::setCalculateNormalVectors,
           "Set whether normal vectors are needed for the supplied velocity "
           "field.")
      .def("setIgnoreVoids", &lsAdvect<T, D>::setIgnoreVoids,
           "Set whether voids in the geometry should be ignored during "
           "advection or not.")
      .def("getAdvectedTime", &lsAdvect<T, D>::getAdvectedTime,
           "Get the time passed during advection.")
      .def("getNumberOfTimeSteps", &lsAdvect<T, D>::getNumberOfTimeSteps,
           "Get how many advection steps were performed after the last apply() "
           "call.")
      .def("getTimeStepRatio", &lsAdvect<T, D>::getTimeStepRatio,
           "Get the time step ratio used for advection.")
      .def("getCalculateNormalVectors",
           &lsAdvect<T, D>::getCalculateNormalVectors,
           "Get whether normal vectors are computed during advection.")
      .def("setIntegrationScheme", &lsAdvect<T, D>::setIntegrationScheme,
           "Set the integration scheme to use during advection.")
      .def("setDissipationAlpha", &lsAdvect<T, D>::setDissipationAlpha,
           "Set the dissipation value to use for Lax Friedrichs integration.")
      // need scoped release since we are calling a python method from
      // parallelised C++ code here
      .def("apply", &lsAdvect<T, D>::apply,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Perform advection.");
  // enums
  pybind11::enum_<lsIntegrationSchemeEnum>(module, "lsIntegrationSchemeEnum")
      .value("ENGQUIST_OSHER_1ST_ORDER",
             lsIntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER)
      .value("ENGQUIST_OSHER_2ND_ORDER",
             lsIntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER)
      .value("LAX_FRIEDRICHS_1ST_ORDER",
             lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER)
      .value("LAX_FRIEDRICHS_2ND_ORDER",
             lsIntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER",
             lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER)
      .value("LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             lsIntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
      .value("LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER",
             lsIntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_2ND_ORDER",
             lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER)
      .value("STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);

  // lsBooleanOperation
  pybind11::class_<lsBooleanOperation<T, D>>(module, "lsBooleanOperation")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsBooleanOperationEnum>())
      .def(pybind11::init<lsDomain<T, D> &, lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsDomain<T, D> &,
                          lsBooleanOperationEnum>())
      // methods
      .def("setLevelset", &lsBooleanOperation<T, D>::setLevelSet,
           "Set levelset on which the boolean operation should be performed.")
      .def("setSecondLevelSet", &lsBooleanOperation<T, D>::setSecondLevelSet,
           "Set second levelset for boolean operation.")
      .def("setBooleanOperation",
           &lsBooleanOperation<T, D>::setBooleanOperation,
           "Set which type of boolean operation should be performed.")
      .def("apply", &lsBooleanOperation<T, D>::apply,
           "Perform the boolean operation.");
  // enums
  pybind11::enum_<lsBooleanOperationEnum>(module, "lsBooleanOperationEnum")
      .value("INTERSECT", lsBooleanOperationEnum::INTERSECT)
      .value("UNION", lsBooleanOperationEnum::UNION)
      .value("RELATIVE_COMPLEMENT", lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .value("INVERT", lsBooleanOperationEnum::INVERT);

  // lsCalculateNormalVectors
  pybind11::class_<lsCalculateNormalVectors<T, D>>(module,
                                                   "lsCalculateNormalVectors")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      // methods
      .def("setLevelSet", &lsCalculateNormalVectors<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &lsCalculateNormalVectors<T, D>::apply,
           "Perform normal vector calculation.");

  // lsCheck
  pybind11::class_<lsCheck<T, D>>(module, "lsCheck")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      // methods
      .def("setLevelSet", &lsCheck<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &lsCheck<T, D>::apply, "Perform check.");

  // lsConvexHull
  pybind11::class_<lsConvexHull<T, D>>(module, "lsConvexHull")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh &, lsPointCloud<T, D> &>())
      // methods
      .def("setMesh", &lsConvexHull<T, D>::setMesh,
           "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &lsConvexHull<T, D>::setPointCloud,
           "Set point cloud used to generate mesh.")
      .def("apply", &lsConvexHull<T, D>::apply, "Perform check.");

  // lsDomain
  pybind11::class_<lsDomain<T, D>>(module, "lsDomain")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<hrleCoordType>())
      .def(pybind11::init<hrleCoordType *, lsDomain<T, D>::BoundaryType *>())
      .def(pybind11::init<hrleCoordType *, lsDomain<T, D>::BoundaryType *,
                          hrleCoordType>())
      .def(pybind11::init<std::vector<hrleCoordType>, std::vector<unsigned>,
                          hrleCoordType>())
      .def(pybind11::init<lsDomain<T, D>::PointValueVectorType, hrleCoordType *,
                          lsDomain<T, D>::BoundaryType *>())
      .def(pybind11::init<lsDomain<T, D>::PointValueVectorType, hrleCoordType *,
                          lsDomain<T, D>::BoundaryType *, hrleCoordType>())
      .def(pybind11::init<const lsDomain<T, D> &>())
      // methods
      .def("deepCopy", &lsDomain<T, D>::deepCopy,
           "Copy lsDomain in this lsDomain.")
      .def("getNumberOfSegments", &lsDomain<T, D>::getNumberOfSegments,
           "Get the number of segments, the level set structure is divided "
           "into.")
      .def("getNumberOfPoints", &lsDomain<T, D>::getNumberOfPoints,
           "Get the number of defined level set values.")
      .def("getLevelSetWidth", &lsDomain<T, D>::getLevelSetWidth,
           "Get the number of layers of level set points around the explicit "
           "surface.")
      .def("setLevelSetWidth", &lsDomain<T, D>::setLevelSetWidth,
           "Set the number of layers of level set points which should be "
           "stored around the explicit surface.")
      .def("clearMetaData", &lsDomain<T, D>::clearMetaData,
           "Clear all metadata stored in the level set.")
      .def("print", &lsDomain<T, D>::print, "Print level set structure.");

  // lsFastAdvect
  pybind11::class_<lsFastAdvect<T, D>>(module, "lsFastAdvect")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &,
                          lsFastAdvectDistribution<hrleCoordType, D> &>())
      // methods
      .def("setLevelSet", &lsFastAdvect<T, D>::setLevelSet,
           "Set levelset to advect.")
      .def(
          "setAdvectionDistribution",
          &lsFastAdvect<T, D>::setAdvectionDistribution,
          "Set advection distribution to use as kernel for the fast advection.")
      .def("apply", &lsFastAdvect<T, D>::apply,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Perform advection.");

  // lsFastAdvectDistributions
  pybind11::class_<lsFastAdvectDistribution<T, D>, PylsFastAdvectDistribution>(
      module, "lsFastAdvectDistribution")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("isInside", &lsFastAdvectDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance",
           &lsFastAdvectDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsFastAdvectDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<lsSphereDistribution<T, D>, lsFastAdvectDistribution<T, D>>(
      module, "lsSphereDistribution")
      // constructors
      .def(pybind11::init<T>())
      // methods
      .def("isInside", &lsSphereDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &lsSphereDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsSphereDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<lsBoxDistribution<T, D>, lsFastAdvectDistribution<T, D>>(
      module, "lsBoxDistribution")
      // constructors
      .def(pybind11::init<const std::array<T, D>>())
      // methods
      .def("isInside", &lsBoxDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &lsBoxDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsBoxDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  // lsExpand
  pybind11::class_<lsExpand<T, D>>(module, "lsExpand")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, int>())
      // methods
      .def("setLevelSet", &lsExpand<T, D>::setLevelSet,
           "Set levelset to expand.")
      .def("setWidth", &lsExpand<T, D>::setWidth, "Set the width to expand to.")
      .def("apply", &lsExpand<T, D>::apply, "Perform expansion.");

  // lsFileFormats
  pybind11::enum_<lsFileFormatEnum>(module, "lsFileFormatEnum")
      .value("VTK_LEGACY", lsFileFormatEnum::VTK_LEGACY)
      .value("VTP", lsFileFormatEnum::VTP)
      .value("VTU", lsFileFormatEnum::VTU);

  // lsFromSurfaceMesh
  pybind11::class_<lsFromSurfaceMesh<T, D>>(module, "lsFromSurfaceMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &, lsMesh &>())
      .def(pybind11::init<lsDomain<T, D> &, lsMesh &, bool>())
      // methods
      .def("setLevelSet", &lsFromSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to read into.")
      .def("setMesh", &lsFromSurfaceMesh<T, D>::setMesh,
           "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles",
           &lsFromSurfaceMesh<T, D>::setRemoveBoundaryTriangles,
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("apply", &lsFromSurfaceMesh<T, D>::apply,
           "Construct a levelset from a surface mesh.");

  // lsFromVolumeMesh
  pybind11::class_<lsFromVolumeMesh<T, D>>(module, "lsFromVolumeMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<std::vector<lsDomain<T, D>> &, lsMesh &>())
      .def(pybind11::init<std::vector<lsDomain<T, D>> &, lsMesh &, bool>())
      // methods
      .def("setLevelSets", &lsFromVolumeMesh<T, D>::setLevelSets,
           "Set levelsets to read into.")
      .def("setMesh", &lsFromVolumeMesh<T, D>::setMesh,
           "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles",
           &lsFromVolumeMesh<T, D>::setRemoveBoundaryTriangles,
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("apply", &lsFromVolumeMesh<T, D>::apply,
           "Construct a levelset from a volume mesh.");

  // lsGeometries
  // lsSphere
  pybind11::class_<lsSphere<T, D>>(module, "lsSphere")
      // constructors
      .def(pybind11::init<const std::vector<T> &, T>())
      // methods
      ;
  // lsPlane
  pybind11::class_<lsPlane<T, D>>(module, "lsPlane")
      // constructors
      .def(pybind11::init<const std::vector<T> &, const std::vector<T> &>())
      // methods
      ;

  // lsBox
  pybind11::class_<lsBox<T, D>>(module, "lsBox")
      // constructors
      .def(pybind11::init<const std::vector<T> &, const std::vector<T> &>())
      // methods
      ;

  // lsPointCloud
  pybind11::class_<lsPointCloud<T, D>>(module, "lsPointCloud")
      // constructors
      .def(pybind11::init<const std::vector<std::vector<T>> &>())
      // methods
      .def("insertNextPoint",
           (void (lsPointCloud<T, D>::*)(const std::vector<T> &)) &
               lsPointCloud<T, D>::insertNextPoint);

  // lsMakeGeometry
  pybind11::class_<lsMakeGeometry<T, D>>(module, "lsMakeGeometry")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsSphere<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsPlane<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, const lsBox<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, lsPointCloud<T, D> &>())
      // methods
      .def("setLevelSet", &lsMakeGeometry<T, D>::setLevelSet,
           "Set the levelset in which to create the geometry.")
      .def("setGeometry", (void (lsMakeGeometry<T, D>::*)(lsSphere<T, D> &)) &
                              lsMakeGeometry<T, D>::setGeometry)

      .def("apply", &lsMakeGeometry<T, D>::apply, "Generate the geometry.");

  // lsPointData
  pybind11::class_<lsPointData>(module, "lsPointData")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("insertNextScalarData",
           (void (lsPointData::*)(const lsPointData::ScalarDataType &,
                                  std::string)) &
               lsPointData::insertNextScalarData,
           pybind11::arg("scalars"), pybind11::arg("label") = "Scalars")
      .def("insertNextVectorData",
           (void (lsPointData::*)(const lsPointData::VectorDataType &,
                                  std::string)) &
               lsPointData::insertNextVectorData,
           pybind11::arg("vectors"), pybind11::arg("label") = "Vectors")
      .def("getScalarDataSize", &lsPointData::getScalarDataSize)
      .def("getVectorDataSize", &lsPointData::getVectorDataSize)
      .def("getScalarData",
           (lsPointData::ScalarDataType * (lsPointData::*)(int)) &
               lsPointData::getScalarData)
      .def("getScalarData",
           (lsPointData::ScalarDataType * (lsPointData::*)(std::string)) &
               lsPointData::getScalarData)
      .def("getScalarDataLabel", &lsPointData::getScalarDataLabel)
      .def("getVectorData",
           (lsPointData::VectorDataType * (lsPointData::*)(int)) &
               lsPointData::getVectorData)
      .def("getVectorData",
           (lsPointData::VectorDataType * (lsPointData::*)(std::string)) &
               lsPointData::getVectorData)
      .def("getVectorDataLabel", &lsPointData::getVectorDataLabel);

  // lsMesh
  pybind11::class_<lsMesh, lsPointData>(module, "lsMesh")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("getNodes",
           (std::vector<std::array<double, 3>> & (lsMesh::*)()) &
               lsMesh::getNodes,
           "Get all nodes of the mesh as a list.")
      .def("getNodes",
           (const std::vector<std::array<double, 3>> &(lsMesh::*)() const) &
               lsMesh::getNodes,
           "Get all nodes of the mesh as a list.")
      .def("getVerticies",
           (std::vector<std::array<unsigned, 1>> & (lsMesh::*)()) &
               lsMesh::getElements<1>,
           "Get a list of verticies of the mesh.")
      .def("getLines",
           (std::vector<std::array<unsigned, 2>> & (lsMesh::*)()) &
               lsMesh::getElements<2>,
           "Get a list of lines of the mesh.")
      .def("getTriangles",
           (std::vector<std::array<unsigned, 3>> & (lsMesh::*)()) &
               lsMesh::getElements<3>,
           "Get a list of verticies of the mesh.")
      .def("getTetras",
           (std::vector<std::array<unsigned, 4>> & (lsMesh::*)()) &
               lsMesh::getElements<4>,
           "Get a list of tetrahedrons of the mesh.")
      .def("getHexas",
           (std::vector<std::array<unsigned, 8>> & (lsMesh::*)()) &
               lsMesh::getElements<8>,
           "Get a list of hexahedrons of the mesh.")
      .def("insertNextNode", &lsMesh::insertNextNode,
           "Insert a node in the mesh.")
      .def("insertNextLine", &lsMesh::insertNextLine,
           "Insert a line in the mesh.")
      .def("insertNextTriangle", &lsMesh::insertNextTriangle,
           "Insert a triangle in the mesh.")
      .def("insertNextTetra", &lsMesh::insertNextTetra,
           "Insert a tetrahedron in the mesh.")
      .def("insertNextHexa", &lsMesh::insertNextHexa,
           "Insert a hexahedron in the mesh.")
      .def("removeDuplicateNodes", &lsMesh::removeDuplicateNodes,
           "Remove nodes which occur twice in the mesh, and replace their IDs "
           "in the mesh elements.")
      .def("print", &lsMesh::print, "Print basic information about the mesh.");

  // lsPrune
  pybind11::class_<lsPrune<T, D>>(module, "lsPrune")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      // methods
      .def("setLevelSet", &lsPrune<T, D>::setLevelSet, "Set levelset to prune.")
      .def("apply", &lsPrune<T, D>::apply, "Perform pruning operation.");

  // lsReduce
  pybind11::class_<lsReduce<T, D>>(module, "lsReduce")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &>())
      .def(pybind11::init<lsDomain<T, D> &, int>())
      .def(pybind11::init<lsDomain<T, D> &, int, bool>())
      // methods
      .def("setLevelSet", &lsReduce<T, D>::setLevelSet,
           "Set levelset to reduce.")
      .def("setWidth", &lsReduce<T, D>::setWidth, "Set the width to reduce to.")
      .def("setNoNewSegment", &lsReduce<T, D>::setNoNewSegment,
           "Set whether the levelset should be segmented anew (balanced across "
           "cores) after reduction.")
      .def("apply", &lsReduce<T, D>::apply, "Perform reduction.");

  // lsToDiskMesh
  pybind11::class_<lsToDiskMesh<T, D>>(module, "lsToDiskMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsDomain<T, D> &, lsMesh &>())
      // methods
      .def("setLevelSet", &lsToDiskMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &lsToDiskMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &lsToDiskMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsToMesh
  pybind11::class_<lsToMesh<T, D>>(module, "lsToMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<const lsDomain<T, D> &, lsMesh &>())
      .def(pybind11::init<const lsDomain<T, D> &, lsMesh &, bool>())
      .def(pybind11::init<const lsDomain<T, D> &, lsMesh &, bool, bool>())
      // methods
      .def("setLevelSet", &lsToMesh<T, D>::setLevelSet, "Set levelset to mesh.")
      .def("setMesh", &lsToMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("setOnlyDefined", &lsToMesh<T, D>::setOnlyDefined,
           "Set whether only defined points should be output to the mesh.")
      .def("setOnlyActive", &lsToMesh<T, D>::setOnlyActive,
           "Set whether only level set points <0.5 should be output.")
      .def("apply", &lsToMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsToSurfaceMesh
  pybind11::class_<lsToSurfaceMesh<T, D>>(module, "lsToSurfaceMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<const lsDomain<T, D> &, lsMesh &>())
      // methods
      .def("setLevelSet", &lsToSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &lsToSurfaceMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("apply", &lsToSurfaceMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsToVoxelMesh
  pybind11::class_<lsToVoxelMesh<T, D>>(module, "lsToVoxelMesh")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh &>())
      .def(pybind11::init<const lsDomain<T, D> &, lsMesh &>())
      .def(pybind11::init<const std::vector<const lsDomain<T, D> *> &,
                          lsMesh &>())
      // methods
      .def("insertNextLevelSet", &lsToVoxelMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the mesh.")
      .def("setMesh", &lsToVoxelMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("apply", &lsToVoxelMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsVelocityField
  pybind11::class_<lsVelocityField<T>, PylsVelocityField>(module,
                                                          "lsVelocityField")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("getScalarVelocity", &lsVelocityField<T>::getScalarVelocity,
           "Return the scalar velocity for a point of material at coordinate "
           "with normal vector normal.")
      .def("getVectorVelocity", &lsVelocityField<T>::getVectorVelocity,
           "Return the vector velocity for a point of material at coordinate "
           "with normal vector normal.")
      .def("getDissipationAlpha", &lsVelocityField<T>::getDissipationAlpha,
           "Return the analytical dissipation alpha value if the "
           "lsLocalLaxFriedrichsAnalytical scheme is used for advection.");

  // lsVTKReader
  pybind11::class_<lsVTKReader>(module, "lsVTKReader")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh &>())
      .def(pybind11::init<lsMesh &, std::string>())
      .def(pybind11::init<lsMesh &, lsFileFormatEnum, std::string>())
      // methods
      .def("setMesh", &lsVTKReader::setMesh, "Set the mesh to read into.")
      .def("setFileFormat", &lsVTKReader::setFileFormat,
           "Set the file format of the file to be read.")
      .def("setFileName", &lsVTKReader::setFileName,
           "Set the name of the input file.")
      .def("apply", &lsVTKReader::apply, "Read the mesh.");

  // lsVTKWriter
  pybind11::class_<lsVTKWriter>(module, "lsVTKWriter")
      // constructors
      .def(pybind11::init<>())
      .def(pybind11::init<lsMesh &>())
      .def(pybind11::init<lsMesh &, std::string>())
      .def(pybind11::init<lsMesh &, lsFileFormatEnum, std::string>())
      // methods
      .def("setMesh", &lsVTKWriter::setMesh, "Set the mesh to output.")
      .def("setFileFormat", &lsVTKWriter::setFileFormat,
           "Set the file format, the mesh should be written to.")
      .def("setFileName", &lsVTKWriter::setFileName,
           "Set the name of the output file.")
      .def("apply", &lsVTKWriter::apply, "Write the mesh.");
}
