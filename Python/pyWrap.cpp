/*
  This file is used to generate the python module of viennals.
  It uses pybind11 to create the modules.

  All necessary headers are included here and the interface
  of the classes which should be exposed defined
*/

// correct module name macro
#define TOKENPASTE_INTERNAL(x, y, z) x##y##z
#define TOKENPASTE(x, y, z) TOKENPASTE_INTERNAL(x, y, z)
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)
#define VIENNALS_MODULE_VERSION STRINGIZE(VIENNALS_VERSION)

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// all header files which define API functions
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateCurvatures.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCheck.hpp>
#include <lsConvexHull.hpp>
#include <lsDetectFeatures.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsExtrude.hpp>
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
#include <lsRemoveStrayPoints.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>
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
                      const vectorType &normalVector,
                      unsigned long pointId) override {
    PYBIND11_OVERLOAD(T, lsVelocityField<T>, getScalarVelocity, coordinate,
                      material, normalVector, pointId);
  }

  vectorType getVectorVelocity(const vectorType &coordinate, int material,
                               const vectorType &normalVector,
                               unsigned long pointId) override {
    PYBIND11_OVERLOAD(vectorType, lsVelocityField<T>, getVectorVelocity,
                      coordinate, material, normalVector, pointId);
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
    PYBIND11_OVERLOAD(bool, ClassType, isInside, initial, candidate, eps);
  }

  T getSignedDistance(const vectorType &initial, const vectorType &candidate,
                      unsigned long initialPointId) const override {
    PYBIND11_OVERLOAD_PURE(T, ClassType, getSignedDistance, initial, candidate,
                           initialPointId);
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
  module.doc() =
      "ViennaLS is a header-only C++ level set library developed for high "
      "performance topography simulations. The main design goals are "
      "simplicity and efficiency, tailored towards scientific simulations. "
      "ViennaLS can also be used for visualisation applications, although this "
      "is not the main design target.";

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

  pybind11::class_<lsCalculateCurvatures<T, D>,
                   lsSmartPointer<lsCalculateCurvatures<T, D>>>(
      module, "lsCalculateCurvatures")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsCalculateCurvatures<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsCalculateCurvatures<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init([](lsSmartPointer<lsDomain<T, D>> &domain,
                             lsCurvatureEnum type) {
        return lsSmartPointer<lsCalculateCurvatures<T, D>>::New(domain, type);
      }))
      // methods
      .def("setLevelSet", &lsCalculateCurvatures<T, D>::setLevelSet,
           "Set levelset for which to calculate the curvatures.")
      .def("setCurvatureType", &lsCalculateCurvatures<T, D>::setCurvatureType,
           "Set which method to use for calculation: Defaults to mean "
           "curvature.")
      .def("setMaxValue", &lsCalculateCurvatures<T, D>::setMaxValue,
           "Curvatures will be calculated for all LS values < maxValue.")
      .def("apply", &lsCalculateCurvatures<T, D>::apply,
           "Perform normal vector calculation.");

  // enums
  pybind11::enum_<lsCurvatureEnum>(module, "lsCurvatureEnum")
      .value("MEAN_CURVATURE", lsCurvatureEnum::MEAN_CURVATURE)
      .value("GAUSSIAN_CURVATURE", lsCurvatureEnum::GAUSSIAN_CURVATURE)
      .value("MEAN_AND_GAUSSIAN_CURVATURE",
             lsCurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE);

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
      .def(pybind11::init(&lsSmartPointer<lsConvexHull<T, D>>::New<
                          lsSmartPointer<lsMesh<T>> &,
                          lsSmartPointer<lsPointCloud<T, D>> &>))
      // methods
      .def("setMesh", &lsConvexHull<T, D>::setMesh,
           "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &lsConvexHull<T, D>::setPointCloud,
           "Set point cloud used to generate mesh.")
      .def("apply", &lsConvexHull<T, D>::apply, "Generate Hull.");

  // lsDetectFeatures
  pybind11::class_<lsDetectFeatures<T, D>,
                   lsSmartPointer<lsDetectFeatures<T, D>>>(module,
                                                           "lsDetectFeatures")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsDetectFeatures<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsDetectFeatures<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsDetectFeatures<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &, T>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init([](lsSmartPointer<lsDomain<T, D>> &domain, T maxValue,
                             lsFeatureDetectionEnum type) {
        return lsSmartPointer<lsDetectFeatures<T, D>>::New(domain, maxValue,
                                                           type);
      }))
      .def("setDetectionThreshold",
           &lsDetectFeatures<T, D>::setDetectionThreshold,
           "Set the curvature value above which a point is considered a "
           "feature.")
      .def("setDetectionMethod", &lsDetectFeatures<T, D>::setDetectionMethod,
           "Set which method to use to detect features. Defaults to Curvature.")
      .def("apply", &lsDetectFeatures<T, D>::apply, "Detect features.");

  // enums
  pybind11::enum_<lsFeatureDetectionEnum>(module, "lsFeatureDetectionEnum")
      .value("CURVATURE", lsFeatureDetectionEnum::CURVATURE)
      .value("NORMALS_ANGLE", lsFeatureDetectionEnum::NORMALS_ANGLE);

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
      // allow filehandle to be passed and default to python standard output
      .def("print", [](lsDomain<T, D>& d, pybind11::object fileHandle) {
          if (!(pybind11::hasattr(fileHandle,"write") &&
          pybind11::hasattr(fileHandle,"flush") )){
               throw pybind11::type_error("MyClass::read_from_file_like_object(file): incompatible function argument:  `file` must be a file-like object, but `"
                                        +(std::string)(pybind11::repr(fileHandle))+"` provided"
               );
          }
          pybind11::detail::pythonbuf buf(fileHandle);
          std::ostream stream(&buf);
          d.print(stream);
          }, pybind11::arg("stream") = pybind11::module::import("sys").attr("stdout"));

  // enums
  pybind11::enum_<lsBoundaryConditionEnum<D>>(module, "lsBoundaryConditionEnum")
      .value("REFLECTIVE_BOUNDARY",
             lsBoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY)
      .value("INFINITE_BOUNDARY", lsBoundaryConditionEnum<D>::INFINITE_BOUNDARY)
      .value("PERIODIC_BOUNDARY", lsBoundaryConditionEnum<D>::PERIODIC_BOUNDARY)
      .value("POS_INFINITE_BOUNDARY",
             lsBoundaryConditionEnum<D>::POS_INFINITE_BOUNDARY)
      .value("NEG_INFINITE_BOUNDARY",
             lsBoundaryConditionEnum<D>::NEG_INFINITE_BOUNDARY);

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

  // lsExtrude
  // Does not work in current implementation, because one can not import both 2D
  // and 3D ViennaLS libraries in Python in the same file
  //   pybind11::class_<lsExtrude<T>, lsSmartPointer<lsExtrude<T>>>(module,
  //                                                                "lsExtrude")
  //       // constructors
  //       .def(pybind11::init(&lsSmartPointer<lsExtrude<T>>::New<>))
  //       .def(
  //           pybind11::init(&lsSmartPointer<lsExtrude<T>>::New<
  //                          lsSmartPointer<lsDomain<T, 2>> &,
  //                          lsSmartPointer<lsDomain<T, 3>> &, std::array<T,
  //                          2>, const int,
  //                          std::array<lsBoundaryConditionEnum<3>, 3>>))
  //       // methods
  //       .def("setInputLevelSet", &lsExtrude<T>::setInputLevelSet,
  //            "Set 2D input Level Set")
  //       .def("setOutputLevelSet", &lsExtrude<T>::setOutputLevelSet,
  //            "Set 3D output Level Set")
  //       .def("setExtent", &lsExtrude<T>::setExtent,
  //            "Set the extent in the extruded dimension")
  //       .def("setExtrudeDimension", &lsExtrude<T>::setExtrudeDimension,
  //            "Set the dimension which should be extruded")
  //       .def("setBoundaryConditions",
  //            pybind11::overload_cast<std::array<lsBoundaryConditionEnum<3>,
  //            3>>(
  //                &lsExtrude<T>::setBoundaryConditions),
  //            "Set the boundary conditions in the 3D extruded domain.")
  //       .def("setBoundaryConditions",
  //            pybind11::overload_cast<lsBoundaryConditionEnum<3> *>(
  //                &lsExtrude<T>::setBoundaryConditions),
  //            "Set the boundary conditions in the 3D extruded domain.")
  //       .def("apply", &lsExtrude<T>::apply, "Perform extrusion.");

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
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(&lsSmartPointer<lsFromSurfaceMesh<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsMesh<T>> &, bool>))
      // methods
      .def("setLevelSet", &lsFromSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to read into.")
      .def("setMesh", &lsFromSurfaceMesh<T, D>::setMesh,
           "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles",
           static_cast<void (lsFromSurfaceMesh<T, D>::*)(bool)>(
               &lsFromSurfaceMesh<T, D>::setRemoveBoundaryTriangles),
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("setRemoveBoundaryTriangles",
           static_cast<void (lsFromSurfaceMesh<T, D>::*)(std::array<bool, 3>)>(
               &lsFromSurfaceMesh<T, D>::setRemoveBoundaryTriangles),
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
                          typename lsFromVolumeMesh<T, D>::GridType &,
                          lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(&lsSmartPointer<lsFromVolumeMesh<T, D>>::New<
                          typename lsFromVolumeMesh<T, D>::GridType &,
                          lsSmartPointer<lsMesh<T>> &, bool>))
      // methods
      .def("setGrid", &lsFromVolumeMesh<T, D>::setGrid,
           "Set the grid used to read in the level sets.")
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
      .def(pybind11::init(&lsSmartPointer<lsSphere<T, D>>::New<
                          const std::vector<T> & /*origin*/, T /*radius*/>),
           pybind11::arg("origin"), pybind11::arg("radius"));
  // lsPlane
  pybind11::class_<lsPlane<T, D>, lsSmartPointer<lsPlane<T, D>>>(module,
                                                                 "lsPlane")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPlane<T, D>>::New<
                          const std::vector<T> & /*origin*/,
                          const std::vector<T> & /*normal*/>),
           pybind11::arg("origin"), pybind11::arg("normal"));
  // lsBox
  pybind11::class_<lsBox<T, D>, lsSmartPointer<lsBox<T, D>>>(module, "lsBox")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsBox<T, D>>::New<
                          const std::vector<T> & /*minPoint*/,
                          const std::vector<T> & /*maxPoint*/>),
           pybind11::arg("minPoint"), pybind11::arg("maxPoint"));
  // lsCylinder
  pybind11::class_<lsCylinder<T, D>, lsSmartPointer<lsCylinder<T, D>>>(
      module, "lsCylinder")
      // constructors
      .def(pybind11::init(
               &lsSmartPointer<lsCylinder<T, D>>::New<
                   const std::vector<T> & /*origin*/,
                   const std::vector<T> & /*axisDirection*/, const T /*height*/,
                   const T /*radius*/, const T /*topRadius*/>),
           pybind11::arg("origin"), pybind11::arg("axisDirection"),
           pybind11::arg("height"), pybind11::arg("radius"),
           pybind11::arg("topRadius") = 0.);
  // lsPointCloud
  pybind11::class_<lsPointCloud<T, D>, lsSmartPointer<lsPointCloud<T, D>>>(
      module, "lsPointCloud")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPointCloud<T, D>>::New<
                          const std::vector<std::vector<T>> &>))
      // methods
      .def("insertNextPoint",
           (void(lsPointCloud<T, D>::*)(const std::vector<T> &)) &
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
                          lsSmartPointer<lsCylinder<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsMakeGeometry<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &,
                          lsSmartPointer<lsPointCloud<T, D>> &>))
      // methods
      .def("setLevelSet", &lsMakeGeometry<T, D>::setLevelSet,
           "Set the levelset in which to create the geometry.")
      .def("setGeometry",
           (void(lsMakeGeometry<T, D>::*)(lsSmartPointer<lsSphere<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(lsMakeGeometry<T, D>::*)(lsSmartPointer<lsPlane<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(lsMakeGeometry<T, D>::*)(lsSmartPointer<lsBox<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(lsMakeGeometry<T, D>::*)(lsSmartPointer<lsCylinder<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(lsMakeGeometry<T, D>::*)(lsSmartPointer<lsPointCloud<T, D>>)) &
               lsMakeGeometry<T, D>::setGeometry)
      .def("setIgnoreBoundaryConditions",
           (void(lsMakeGeometry<T, D>::*)(bool)) &
               lsMakeGeometry<T, D>::setIgnoreBoundaryConditions)
      .def("setIgnoreBoundaryConditions",
           (void(lsMakeGeometry<T, D>::*)(std::array<bool, 3>)) &
               lsMakeGeometry<T, D>::setIgnoreBoundaryConditions)
      .def("apply", &lsMakeGeometry<T, D>::apply, "Generate the geometry.");

  // lsMarkVoidPoints
  pybind11::class_<lsMarkVoidPoints<T, D>,
                   lsSmartPointer<lsMarkVoidPoints<T, D>>>(module,
                                                           "lsMarkVoidPoints")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsMarkVoidPoints<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsMarkVoidPoints<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      .def(pybind11::init(&lsSmartPointer<lsMarkVoidPoints<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &, bool &>))
      // methods
      .def("setLevelSet", &lsMarkVoidPoints<T, D>::setLevelSet,
           "Set the levelset to mark void points in.")
      .def("setReverseVoidDetection",
           &lsMarkVoidPoints<T, D>::setReverseVoidDetection,
           "Reverse the logic of detecting the top surface.")
      .def("setDetectLargestSurface",
           &lsMarkVoidPoints<T, D>::setDetectLargestSurface,
           "Set that the top surface should be the one with the most connected "
           "LS points.")
      .def("setVoidTopSurface", &lsMarkVoidPoints<T, D>::setVoidTopSurface,
           "Set the logic by which to choose the surface which is non-void. "
           "All other connected surfaces will then be marked as void points.")
      .def("setSaveComponentsId", &lsMarkVoidPoints<T, D>::setSaveComponentIds,
           "Save the connectivity information of all LS points in the "
           "pointData of the level set.")
      .def("apply", &lsMarkVoidPoints<T, D>::apply, "Mark void points.");

  // lsVoidTopSurfaceEnum
  pybind11::enum_<lsVoidTopSurfaceEnum>(module, "lsVoidTopSurfaceEnum")
      .value("LEX_LOWEST", lsVoidTopSurfaceEnum::LEX_LOWEST)
      .value("LEX_HIGHEST", lsVoidTopSurfaceEnum::LEX_HIGHEST)
      .value("LARGEST", lsVoidTopSurfaceEnum::LARGEST)
      .value("SMALLEST", lsVoidTopSurfaceEnum::SMALLEST);

  // lsPointData
  pybind11::class_<lsPointData<T>, lsSmartPointer<lsPointData<T>>>(
      module, "lsPointData")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsPointData<T>>::New<>))
      // methods
      .def("insertNextScalarData",
           (void(lsPointData<T>::*)(const lsPointData<T>::ScalarDataType &,
                                    std::string)) &
               lsPointData<T>::insertNextScalarData,
           pybind11::arg("scalars"), pybind11::arg("label") = "Scalars")
      .def("insertNextVectorData",
           (void(lsPointData<T>::*)(const lsPointData<T>::VectorDataType &,
                                    std::string)) &
               lsPointData<T>::insertNextVectorData,
           pybind11::arg("vectors"), pybind11::arg("label") = "Vectors")
      .def("getScalarDataSize", &lsPointData<T>::getScalarDataSize)
      .def("getVectorDataSize", &lsPointData<T>::getVectorDataSize)
      .def("getScalarData",
           (lsPointData<T>::ScalarDataType * (lsPointData<T>::*)(int)) &
               lsPointData<T>::getScalarData)
      .def("getScalarData",
           (lsPointData<T>::ScalarDataType * (lsPointData<T>::*)(std::string)) &
               lsPointData<T>::getScalarData)
      .def("getScalarDataLabel", &lsPointData<T>::getScalarDataLabel)
      .def("getVectorData",
           (lsPointData<T>::VectorDataType * (lsPointData<T>::*)(int)) &
               lsPointData<T>::getVectorData)
      .def("getVectorData",
           (lsPointData<T>::VectorDataType * (lsPointData<T>::*)(std::string)) &
               lsPointData<T>::getVectorData)
      .def("getVectorDataLabel", &lsPointData<T>::getVectorDataLabel);

  // lsMesh<T>
  pybind11::class_<lsMesh<T>, lsSmartPointer<lsMesh<T>>>(module, "lsMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsMesh<T>>::New<>))
      // methods
      .def("getNodes",
           (std::vector<std::array<double, 3>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getNodes,
           "Get all nodes of the mesh as a list.")
      .def("getNodes",
           (std::vector<std::array<double, 3>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getNodes,
           "Get all nodes of the mesh as a list.")
      .def("getVerticies",
           (std::vector<std::array<unsigned, 1>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getElements<1>,
           "Get a list of verticies of the mesh.")
      .def("getLines",
           (std::vector<std::array<unsigned, 2>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getElements<2>,
           "Get a list of lines of the mesh.")
      .def("getTriangles",
           (std::vector<std::array<unsigned, 3>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getElements<3>,
           "Get a list of triangles of the mesh.")
      .def("getTetras",
           (std::vector<std::array<unsigned, 4>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getElements<4>,
           "Get a list of tetrahedrons of the mesh.")
      .def("getHexas",
           (std::vector<std::array<unsigned, 8>> & (lsMesh<T>::*)()) &
               lsMesh<T>::getElements<8>,
           "Get a list of hexahedrons of the mesh.")
      .def("getPointData",
           (lsPointData<T> & (lsMesh<T>::*)()) & lsMesh<T>::getPointData,
           "Return a reference to the point data of the mesh.")
      .def("getCellData",
           (lsPointData<T> & (lsMesh<T>::*)()) & lsMesh<T>::getCellData,
           "Return a reference to the cell data of the mesh.")
      .def("insertNextNode", &lsMesh<T>::insertNextNode,
           "Insert a node in the mesh.")
      .def("insertNextVertex", &lsMesh<T>::insertNextVertex,
           "Insert a vertex in the mesh.")
      .def("insertNextLine", &lsMesh<T>::insertNextLine,
           "Insert a line in the mesh.")
      .def("insertNextTriangle", &lsMesh<T>::insertNextTriangle,
           "Insert a triangle in the mesh.")
      .def("insertNextTetra", &lsMesh<T>::insertNextTetra,
           "Insert a tetrahedron in the mesh.")
      .def("insertNextHexa", &lsMesh<T>::insertNextHexa,
           "Insert a hexahedron in the mesh.")
      .def("removeDuplicateNodes", &lsMesh<T>::removeDuplicateNodes,
           "Remove nodes which occur twice in the mesh, and replace their IDs "
           "in the mesh elements.")
      .def("print", &lsMesh<T>::print,
           "Print basic information about the mesh.");

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

  // lsRemoveStrayPoints
  pybind11::class_<lsRemoveStrayPoints<T, D>,
                   lsSmartPointer<lsRemoveStrayPoints<T, D>>>(
      module, "lsRemoveStrayPoints")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsRemoveStrayPoints<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsRemoveStrayPoints<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // methods
      .def("setLevelSet", &lsRemoveStrayPoints<T, D>::setLevelSet,
           "Set levelset for stray point removal.")
      .def("setVoidTopSurface", &lsRemoveStrayPoints<T, D>::setVoidTopSurface,
           "Set the logic by which to choose the surface which should be kept. "
           "All other LS values will be marked as stray points and removed.")
      .def("apply", &lsRemoveStrayPoints<T, D>::apply, "Remove stray points.");

  // lsToDiskMesh
  pybind11::class_<lsToDiskMesh<T, D>, lsSmartPointer<lsToDiskMesh<T, D>>>(
      module, "lsToDiskMesh")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsToDiskMesh<T, D>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsToDiskMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh<T>> &>))
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
                                               lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsToMesh<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               lsSmartPointer<lsMesh<T>> &,
                                               bool>))
      .def(pybind11::init(
          &lsSmartPointer<lsToMesh<T, D>>::New<lsSmartPointer<lsDomain<T, D>> &,
                                               lsSmartPointer<lsMesh<T>> &,
                                               bool, bool>))
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
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh<T>> &>))
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
      .def(pybind11::init(&lsSmartPointer<lsToVoxelMesh<T, D>>::New<
                          lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsToVoxelMesh<T, D>>::New<
              lsSmartPointer<lsDomain<T, D>> &, lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(&lsSmartPointer<lsToVoxelMesh<T, D>>::New<
                          std::vector<lsSmartPointer<lsDomain<T, D>>> &,
                          lsSmartPointer<lsMesh<T>> &>))
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
  pybind11::class_<lsVTKReader<T>, lsSmartPointer<lsVTKReader<T>>>(
      module, "lsVTKReader")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsVTKReader<T>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKReader<T>>::New<lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKReader<T>>::New<lsSmartPointer<lsMesh<T>> &,
                                               std::string>))
      .def(pybind11::init([](lsSmartPointer<lsMesh<T>> &mesh,
                             lsFileFormatEnum format, std::string s) {
        return lsSmartPointer<lsVTKReader<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &lsVTKReader<T>::setMesh, "Set the mesh to read into.")
      .def("setFileFormat", &lsVTKReader<T>::setFileFormat,
           "Set the file format of the file to be read.")
      .def("setFileName", &lsVTKReader<T>::setFileName,
           "Set the name of the input file.")
      .def("apply", &lsVTKReader<T>::apply, "Read the mesh.");

  // lsVTKWriter
  pybind11::class_<lsVTKWriter<T>, lsSmartPointer<lsVTKWriter<T>>>(
      module, "lsVTKWriter")
      // constructors
      .def(pybind11::init(&lsSmartPointer<lsVTKWriter<T>>::New<>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKWriter<T>>::New<lsSmartPointer<lsMesh<T>> &>))
      .def(pybind11::init(
          &lsSmartPointer<lsVTKWriter<T>>::New<lsSmartPointer<lsMesh<T>> &,
                                               std::string>))
      .def(pybind11::init([](lsSmartPointer<lsMesh<T>> &mesh,
                             lsFileFormatEnum format, std::string s) {
        return lsSmartPointer<lsVTKWriter<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &lsVTKWriter<T>::setMesh, "Set the mesh to output.")
      .def("setFileFormat", &lsVTKWriter<T>::setFileFormat,
           "Set the file format, the mesh should be written to.")
      .def("setFileName", &lsVTKWriter<T>::setFileName,
           "Set the name of the output file.")
      .def("apply", &lsVTKWriter<T>::apply, "Write the mesh.");

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

// lsWriteVisualizationMesh
#ifdef VIENNALS_USE_VTK
  pybind11::class_<lsWriteVisualizationMesh<T, D>,
                   lsSmartPointer<lsWriteVisualizationMesh<T, D>>>(
      module, "lsWriteVisualizationMesh")
      // constructors
      .def(pybind11::init(
          &lsSmartPointer<lsWriteVisualizationMesh<T, D>>::New<>))
      .def(pybind11::init(&lsSmartPointer<lsWriteVisualizationMesh<T, D>>::New<
                          lsSmartPointer<lsDomain<T, D>> &>))
      // methods
      .def("insertNextLevelSet",
           &lsWriteVisualizationMesh<T, D>::insertNextLevelSet,
           "Insert next level set to convert. Bigger level sets wrapping "
           "smaller ones, should be inserted last.")
      .def("setFileName", &lsWriteVisualizationMesh<T, D>::setFileName,
           "Set Name of File to write.")
      .def("setExtractHullMesh",
           &lsWriteVisualizationMesh<T, D>::setExtractHullMesh,
           "Whether to extract a hull mesh. Defaults to false.")
      .def("setExtractVolumeMesh",
           &lsWriteVisualizationMesh<T, D>::setExtractVolumeMesh,
           " Whether to extract a tetra volume mesh. Defaults to true.")
      .def("apply", &lsWriteVisualizationMesh<T, D>::apply,
           "Make and write mesh.");
#endif
}
