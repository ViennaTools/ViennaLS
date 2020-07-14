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
#include <lsFileFormats.hpp>
#include <lsFromMesh.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsMesh.hpp>
#include <lsPointData.hpp>
#include <lsPrune.hpp>
#include <lsReader.hpp>
#include <lsReduce.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriter.hpp>

// always use double for python export
typedef double T;
// get dimension from cmake define
constexpr int D = VIENNALS_PYTHON_DIMENSION;

PYBIND11_DECLARE_HOLDER_TYPE(TemplateType, lsSmartPointer<TemplateType>);

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

// lsGeometricAdvectDistribution
class PylsGeometricAdvectDistribution
    : public lsGeometricAdvectDistribution<T, D> {
  typedef std::array<hrleCoordType, 3> vectorType;
  typedef std::array<hrleCoordType, 6> boundsType;
  typedef lsGeometricAdvectDistribution<T, D> ClassType;
  using lsGeometricAdvectDistribution<T, D>::lsGeometricAdvectDistribution;

public:
  bool isInside(const vectorType &initial, const vectorType &candidate,
                double eps = 0.) const override {
    PYBIND11_OVERLOAD_PURE(bool, ClassType, isInside, initial, candidate, eps);
  }

  T getSignedDistance(const vectorType &initial,
                      const vectorType &candidate) const override {
    PYBIND11_OVERLOAD_PURE(T, ClassType, getSignedDistance, initial, candidate);
  }

  boundsType getBounds() const override {
    PYBIND11_OVERLOAD_PURE(boundsType, ClassType, getBounds);
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

  // wrap omp_set_num_threads to control number of threads
  module.def("setNumThreads", &omp_set_num_threads);

  // lsAdvect
  pybind11::class_<lsAdvect<T, D>, lsSmartPointer<lsAdvect<T, D>>>(module,
                                                                   "lsAdvect")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsAdvect<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsAdvect<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsAdvect<T, D>>::New<
                          lsSmartPointer<lsVelocityField<T>> &>))
      .def(pybind11::init(&lsSmartPointer<lsAdvect<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsVelocityField<T>> &>))
      // getters and setters
      .def("insertNextLevelSet", &lsAdvect<T, D>::insertNextLevelSet,
           "Insert next level set to use for advection.")
      .def("setVelocityField",
           &lsAdvect<T, D>::setVelocityField<PylsVelocityField>,
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
  pybind11::class_<lsBooleanOperation<T, D>,
                   lsSmartPointer<lsBooleanOperation<T, D>>>(
      module, "lsBooleanOperation")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsBooleanOperation<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsBooleanOperation<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsBooleanOperation<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsDomain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init([](lsSmartPointer<lsDomain<T, D>> &domain,
                             lsBooleanOperationEnum op) {
        return lsSmartPointer<lsBooleanOperation<T, D>>::New(domain, op);
      }))
      .def(pybind11::init([](lsSmartPointer<lsDomain<T, D>> &domainA,
                             lsSmartPointer<lsDomain<T, D>> &domainB,
                             lsBooleanOperationEnum op) {
        return lsSmartPointer<lsBooleanOperation<T, D>>::New(domainA, domainB,
                                                             op);
      }))
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
  pybind11::class_<lsCalculateNormalVectors<T, D>,
                   lsSmartPointer<lsCalculateNormalVectors<T, D>>>(
      module, "lsCalculateNormalVectors")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsCalculateNormalVectors<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsCalculateNormalVectors<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // methods
      .def("setLevelSet", &lsCalculateNormalVectors<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &lsCalculateNormalVectors<T, D>::apply,
           "Perform normal vector calculation.");

  // lsCheck
  pybind11::class_<lsCheck<T, D>, lsSmartPointer<lsCheck<T, D>>>(module,
                                                                 "lsCheck")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsCheck<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsCheck<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // methods
      .def("setLevelSet", &lsCheck<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &lsCheck<T, D>::apply, "Perform check.");

  // lsConvexHull
  pybind11::class_<lsConvexHull<T, D>, lsSmartPointer<lsConvexHull<T, D>>>(
      module, "lsConvexHull")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsConvexHull<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsConvexHull<T, D>>::New<
              lsSmartPointer<lsMesh> &, lsSmartPointer<lsPointCloud<T, D>> &>))
      // methods
      .def("setMesh", &lsConvexHull<T, D>::setMesh,
           "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &lsConvexHull<T, D>::setPointCloud,
           "Set point cloud used to generate mesh.")
      .def("apply", &lsConvexHull<T, D>::apply, "Perform check.");

  // lsDomain
  pybind11::class_<lsDomain<T, D>, lsSmartPointer<lsDomain<T, D>>>(module,
                                                                   "lsDomain")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsDomain<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsDomain<T, D>>::New<hrleCoordType>))
      .def(pybind11::init(
          &lsSmartPointer<lsDomain<T, D>>::New<hrleCoordType *,
                                               lsDomain<T, D>::BoundaryType *>))
      .def(pybind11::init(
          &lsSmartPointer<lsDomain<T, D>>::New<
              hrleCoordType *, lsDomain<T, D>::BoundaryType *, hrleCoordType>))
      .def(pybind11::init(
          &lsSmartPointer<lsDomain<T, D>>::New<std::vector<hrleCoordType>,
                                               std::vector<unsigned>,
                                               hrleCoordType>))
      .def(pybind11::init(&lsSmartPointer<lsDomain<T, D>>::New<
                          lsDomain<T, D>::PointValueVectorType, hrleCoordType *,
                          lsDomain<T, D>::BoundaryType *>))
      .def(pybind11::init(&lsSmartPointer<lsDomain<T, D>>::New<
                          lsDomain<T, D>::PointValueVectorType, hrleCoordType *,
                          lsDomain<T, D>::BoundaryType *, hrleCoordType>))
      .def(pybind11::init(&lsSmartPointer<lsDomain<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
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

  // lsGeometricAdvect
  pybind11::class_<lsGeometricAdvect<T, D>,
                   lsSmartPointer<lsGeometricAdvect<T, D>>>(module,
                                                            "lsGeometricAdvect")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsGeometricAdvect<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsGeometricAdvect<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &,
              lsSmartPointer<lsGeometricAdvectDistribution<hrleCoordType, D>>
                  &>))
      // methods
      .def("setLevelSet", &lsGeometricAdvect<T, D>::setLevelSet,
           "Set levelset to advect.")
      .def(
          "setAdvectionDistribution",
          &lsGeometricAdvect<T, D>::setAdvectionDistribution<
              PylsGeometricAdvectDistribution>,
          "Set advection distribution to use as kernel for the fast advection.")
      .def("apply", &lsGeometricAdvect<T, D>::apply,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Perform advection.");

  // lsGeometricAdvectDistributions
  pybind11::class_<lsGeometricAdvectDistribution<T, D>,
                   lsSmartPointer<lsGeometricAdvectDistribution<T, D>>,
                   PylsGeometricAdvectDistribution>(
      module, "lsGeometricAdvectDistribution")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("isInside", &lsGeometricAdvectDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance",
           &lsGeometricAdvectDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsGeometricAdvectDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<lsSphereDistribution<T, D>,
                   lsSmartPointer<lsSphereDistribution<T, D>>,
                   lsGeometricAdvectDistribution<T, D>>(module,
                                                        "lsSphereDistribution")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsSphereDistribution<T, D>>::New<T, T>))
      // methods
      .def("isInside", &lsSphereDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &lsSphereDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsSphereDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<lsBoxDistribution<T, D>,
                   lsSmartPointer<lsBoxDistribution<T, D>>,
                   lsGeometricAdvectDistribution<T, D>>(module,
                                                        "lsBoxDistribution")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsBoxDistribution<T, D>>::New<const std::array<T, 3>,
                                                        T>))
      // methods
      .def("isInside", &lsBoxDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &lsBoxDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &lsBoxDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  // lsExpand
  pybind11::class_<lsExpand<T, D>, lsSmartPointer<lsExpand<T, D>>>(module,
                                                                   "lsExpand")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsExpand<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsExpand<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsExpand<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               int>))
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
  pybind11::class_<lsFromSurfaceMesh<T, D>,
                   lsSmartPointer<lsFromSurfaceMesh<T, D>>>(module,
                                                            "lsFromSurfaceMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsFromSurfaceMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsFromSurfaceMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(&lsSmartPointer<lsFromSurfaceMesh<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsMesh> &, bool>))
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
  pybind11::class_<lsFromVolumeMesh<T, D>,
                   lsSmartPointer<lsFromVolumeMesh<T, D>>>(module,
                                                           "lsFromVolumeMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsFromVolumeMesh<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsFromVolumeMesh<T, D>>::New<
                          std::vector<lsSmartPointer<lsDomain<T, D>>> &,
                          lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(&lsSmartPointer<lsFromVolumeMesh<T, D>>::New<
                          std::vector<lsSmartPointer<lsDomain<T, D>>> &,
                          lsSmartPointer<lsMesh> &, bool>))
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
  pybind11::class_<lsSphere<T, D>, lsSmartPointer<lsSphere<T, D>>>(module,
                                                                   "lsSphere")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsSphere<T, D>>::New<const std::vector<T> &, T>))
      // methods
      ;
  // lsPlane
  pybind11::class_<lsPlane<T, D>, lsSmartPointer<lsPlane<T, D>>>(module,
                                                                 "lsPlane")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsPlane<T, D>>::New<const std::vector<T> &,
                                              const std::vector<T> &>))
      // methods
      ;

  // lsBox
  pybind11::class_<lsBox<T, D>, lsSmartPointer<lsBox<T, D>>>(module, "lsBox")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsBox<T, D>>::New<const std::vector<T> &,
                                            const std::vector<T> &>))
      // methods
      ;

  // lsPointCloud
  pybind11::class_<lsPointCloud<T, D>, lsSmartPointer<lsPointCloud<T, D>>>(
      module, "lsPointCloud")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPointCloud<T, D>>::New<
                          const std::vector<std::vector<T>> &>))
      // methods
      .def("insertNextPoint",
           (void (lsPointCloud<T, D>::*)(const std::vector<T> &)) &
               lsPointCloud<T, D>::insertNextPoint);

  // lsMakeGeometry
  pybind11::class_<lsMakeGeometry<T, D>, lsSmartPointer<lsMakeGeometry<T, D>>>(
      module, "lsMakeGeometry")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsSphere<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsPlane<T, D>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsMakeGeometry<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsBox<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsPointCloud<T, D>> &>))
      // methods
      .def("setLevelSet", &lsMakeGeometry<T, D>::setLevelSet,
           "Set the levelset in which to create the geometry.")
      .def("setGeometry",
           (void (lsMakeGeometry<T, D>::*)(lsSmartPointer<lsSphere<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void (lsMakeGeometry<T, D>::*)(lsSmartPointer<lsPlane<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void (lsMakeGeometry<T, D>::*)(lsSmartPointer<lsBox<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry", (void (lsMakeGeometry<T, D>::*)(
                              lsSmartPointer<lsPointCloud<T, D>>)) &
                              lsMakeGeometry<T, D>::setGeometry)
      .def("apply", &lsMakeGeometry<T, D>::apply, "Generate the geometry.");

  // lsPointData
  pybind11::class_<lsPointData, lsSmartPointer<lsPointData>>(module,
                                                             "lsPointData")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPointData>::New<>))
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
  pybind11::class_<lsMesh, lsSmartPointer<lsMesh>, lsPointData>(module,
                                                                "lsMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsMesh>::New<>))
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
  pybind11::class_<lsPrune<T, D>, lsSmartPointer<lsPrune<T, D>>>(module,
                                                                 "lsPrune")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPrune<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsPrune<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // methods
      .def("setLevelSet", &lsPrune<T, D>::setLevelSet, "Set levelset to prune.")
      .def("apply", &lsPrune<T, D>::apply, "Perform pruning operation.");

  // lsReader
  pybind11::class_<lsReader<T, D>, lsSmartPointer<lsReader<T, D>>>(module,
                                                                   "lsReader")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsReader<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsReader<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsReader<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               std::string>))
      // methods
      .def("setLevelSet", &lsReader<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &lsReader<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &lsReader<T, D>::apply, "Write to file.");

  // lsReduce
  pybind11::class_<lsReduce<T, D>, lsSmartPointer<lsReduce<T, D>>>(module,
                                                                   "lsReduce")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsReduce<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsReduce<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsReduce<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               int>))
      .def(pybind11::init(
          &lsSmartPointer<lsReduce<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               int, bool>))
      // methods
      .def("setLevelSet", &lsReduce<T, D>::setLevelSet,
           "Set levelset to reduce.")
      .def("setWidth", &lsReduce<T, D>::setWidth, "Set the width to reduce to.")
      .def("setNoNewSegment", &lsReduce<T, D>::setNoNewSegment,
           "Set whether the levelset should be segmented anew (balanced across "
           "cores) after reduction.")
      .def("apply", &lsReduce<T, D>::apply, "Perform reduction.");

  // lsToDiskMesh
  pybind11::class_<lsToDiskMesh<T, D>, lsSmartPointer<lsToDiskMesh<T, D>>>(
      module, "lsToDiskMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsToDiskMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsToDiskMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh> &>))
      // methods
      .def("setLevelSet", &lsToDiskMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &lsToDiskMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &lsToDiskMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsToMesh
  pybind11::class_<lsToMesh<T, D>, lsSmartPointer<lsToMesh<T, D>>>(module,
                                                                   "lsToMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsToMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsToMesh<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsToMesh<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               lsSmartPointer<lsMesh> &, bool>))
      .def(pybind11::init(
          &lsSmartPointer<lsToMesh<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               lsSmartPointer<lsMesh> &, bool,
                                               bool>))
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
  pybind11::class_<lsToSurfaceMesh<T, D>,
                   lsSmartPointer<lsToSurfaceMesh<T, D>>>(module,
                                                          "lsToSurfaceMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsToSurfaceMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsToSurfaceMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh> &>))
      // methods
      .def("setLevelSet", &lsToSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &lsToSurfaceMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("apply", &lsToSurfaceMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsToVoxelMesh
  pybind11::class_<lsToVoxelMesh<T, D>, lsSmartPointer<lsToVoxelMesh<T, D>>>(
      module, "lsToVoxelMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsToVoxelMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsToVoxelMesh<T, D>>::New<lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsToVoxelMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(&lsSmartPointer<lsToVoxelMesh<T, D>>::New<
                          std::vector<lsSmartPointer<lsDomain<T, D>>> &,
                          lsSmartPointer<lsMesh> &>))
      // methods
      .def("insertNextLevelSet", &lsToVoxelMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the mesh.")
      .def("setMesh", &lsToVoxelMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("apply", &lsToVoxelMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // lsVelocityField
  pybind11::class_<lsVelocityField<T>, lsSmartPointer<lsVelocityField<T>>,
                   PylsVelocityField>(module, "lsVelocityField")
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
  pybind11::class_<lsVTKReader, lsSmartPointer<lsVTKReader>>(module,
                                                             "lsVTKReader")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsVTKReader>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKReader>::New<lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKReader>::New<lsSmartPointer<lsMesh> &,
                                            std::string>))
      .def(pybind11::init([](lsSmartPointer<lsMesh> &mesh,
                             lsFileFormatEnum format, std::string s) {
        return lsSmartPointer<lsVTKReader>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &lsVTKReader::setMesh, "Set the mesh to read into.")
      .def("setFileFormat", &lsVTKReader::setFileFormat,
           "Set the file format of the file to be read.")
      .def("setFileName", &lsVTKReader::setFileName,
           "Set the name of the input file.")
      .def("apply", &lsVTKReader::apply, "Read the mesh.");

  // lsVTKWriter
  pybind11::class_<lsVTKWriter, lsSmartPointer<lsVTKWriter>>(module,
                                                             "lsVTKWriter")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsVTKWriter>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKWriter>::New<lsSmartPointer<lsMesh> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKWriter>::New<lsSmartPointer<lsMesh> &,
                                            std::string>))
      .def(pybind11::init([](lsSmartPointer<lsMesh> &mesh,
                             lsFileFormatEnum format, std::string s) {
        return lsSmartPointer<lsVTKWriter>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &lsVTKWriter::setMesh, "Set the mesh to output.")
      .def("setFileFormat", &lsVTKWriter::setFileFormat,
           "Set the file format, the mesh should be written to.")
      .def("setFileName", &lsVTKWriter::setFileName,
           "Set the name of the output file.")
      .def("apply", &lsVTKWriter::apply, "Write the mesh.");

  // lsWriter
  pybind11::class_<lsWriter<T, D>, lsSmartPointer<lsWriter<T, D>>>(module,
                                                                   "lsWriter")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsWriter<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsWriter<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsWriter<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               std::string>))
      // methods
      .def("setLevelSet", &lsWriter<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &lsWriter<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &lsWriter<T, D>::apply, "Write to file.");
}
