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
#include <lsCalculateVisibilities.hpp>
#include <lsCheck.hpp>
#include <lsCompareArea.hpp>
#include <lsCompareChamfer.hpp>
#include <lsCompareCriticalDimensions.hpp>
#include <lsCompareNarrowBand.hpp>
#include <lsCompareSparseField.hpp>
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
#include <lsMaterialMap.hpp>
#include <lsMesh.hpp>
#include <lsPointData.hpp>
#include <lsPrune.hpp>
#include <lsReader.hpp>
#include <lsReduce.hpp>
#include <lsRemoveStrayPoints.hpp>
#include <lsSlice.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsTransformMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <lsWriter.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

using namespace viennals;

// always use double for python export
typedef double T;
// get dimension from cmake define
constexpr int D = VIENNALS_PYTHON_DIMENSION;

PYBIND11_DECLARE_HOLDER_TYPE(TemplateType, SmartPointer<TemplateType>);

// define trampoline classes for interface functions
// ALSO NEED TO ADD TRAMPOLINE CLASSES FOR CLASSES
// WHICH HOLD REFERENCES TO INTERFACE(ABSTRACT) CLASSES

// BASE CLASS WRAPPERS
// lsVelocityField only defines interface and has no functionality
class PylsVelocityField : public VelocityField<T> {
  typedef std::array<T, 3> vectorType;
  using VelocityField<T>::VelocityField;

public:
  T getScalarVelocity(const vectorType &coordinate, int material,
                      const vectorType &normalVector,
                      unsigned long pointId) override {
    PYBIND11_OVERLOAD(T, VelocityField<T>, getScalarVelocity, coordinate,
                      material, normalVector, pointId);
  }

  vectorType getVectorVelocity(const vectorType &coordinate, int material,
                               const vectorType &normalVector,
                               unsigned long pointId) override {
    PYBIND11_OVERLOAD(vectorType, VelocityField<T>, getVectorVelocity,
                      coordinate, material, normalVector, pointId);
  }
};

// GeometricAdvectDistribution
class PylsGeometricAdvectDistribution
    : public GeometricAdvectDistribution<T, D> {
  typedef std::array<viennahrle::CoordType, 3> vectorType;
  typedef std::array<viennahrle::CoordType, 6> boundsType;
  typedef GeometricAdvectDistribution<T, D> ClassType;
  using GeometricAdvectDistribution<T, D>::GeometricAdvectDistribution;

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
// class PylsAdvect : public Advect<T, D> {
//   pybind11::object pyObj;
// public:
//   PylsAdvect(Domain<T, D> &passedDomain, VelocityField<T>
//   &passedVelocities) : Advect<T,D>(passedDomain, passedVelocities),
//   pyObj(pybind11::cast(passedVelocities)) {}
//
//   PylsAdvect(VelocityField<T> &passedVelocities) :
//   Advect<T,D>(passedVelocities), pyObj(pybind11::cast(passedVelocities)) {}
// };

// module specification
PYBIND11_MODULE(VIENNALS_MODULE_NAME, module) {
  module.doc() =
      "ViennaLS is a header-only C++ level set library developed for high "
      "performance topography simulations. The main design goals are "
      "simplicity and efficiency, tailored towards scientific simulations. "
      "ViennaLS can also be used for visualization applications, although this "
      "is not the main design target.";

  // set version string of python module
  module.attr("__version__") =
      VIENNALS_MODULE_VERSION; // for some reason this string does not show
  module.attr("version") = VIENNALS_MODULE_VERSION;

  // wrap omp_set_num_threads to control number of threads
  module.def("setNumThreads", &omp_set_num_threads);

  // Logger
  pybind11::enum_<LogLevel>(module, "LogLevel", pybind11::module_local())
      .value("ERROR", LogLevel::ERROR)
      .value("WARNING", LogLevel::WARNING)
      .value("INFO", LogLevel::INFO)
      .value("TIMING", LogLevel::TIMING)
      .value("INTERMEDIATE", LogLevel::INTERMEDIATE)
      .value("DEBUG", LogLevel::DEBUG);

  pybind11::class_<Logger, SmartPointer<Logger>>(module, "Logger",
                                                 pybind11::module_local())
      .def_static("setLogLevel", &Logger::setLogLevel)
      .def_static("getLogLevel", &Logger::getLogLevel)
      .def_static("getInstance", &Logger::getInstance,
                  pybind11::return_value_policy::reference)
      .def("addDebug", &Logger::addDebug)
      .def("addTiming", (Logger & (Logger::*)(const std::string &, double)) &
                            Logger::addTiming)
      .def("addTiming",
           (Logger & (Logger::*)(const std::string &, double, double)) &
               Logger::addTiming)
      .def("addInfo", &Logger::addInfo)
      .def("addWarning", &Logger::addWarning)
      .def("addError", &Logger::addError, pybind11::arg("s"),
           pybind11::arg("shouldAbort") = true)
      .def("print", [](Logger &instance) { instance.print(std::cout); });

  // Advect
  pybind11::class_<Advect<T, D>, SmartPointer<Advect<T, D>>>(module, "Advect")
      // constructors
      .def(pybind11::init(&SmartPointer<Advect<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Advect<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Advect<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           SmartPointer<VelocityField<T>> &>))
      // getters and setters
      .def("insertNextLevelSet", &Advect<T, D>::insertNextLevelSet,
           "Insert next level set to use for advection.")
      .def("setVelocityField", &Advect<T, D>::setVelocityField,
           "Set the velocity to use for advection.")
      .def("setAdvectionTime", &Advect<T, D>::setAdvectionTime,
           "Set the time until when the level set should be advected.")
      .def("setTimeStepRatio", &Advect<T, D>::setTimeStepRatio,
           "Set the maximum time step size relative to grid size. Advection is "
           "only stable for <0.5.")
      .def("setCalculateNormalVectors",
           &Advect<T, D>::setCalculateNormalVectors,
           "Set whether normal vectors are needed for the supplied velocity "
           "field.")
      .def("setIgnoreVoids", &Advect<T, D>::setIgnoreVoids,
           "Set whether voids in the geometry should be ignored during "
           "advection or not.")
      .def("getAdvectedTime", &Advect<T, D>::getAdvectedTime,
           "Get the time passed during advection.")
      .def("getNumberOfTimeSteps", &Advect<T, D>::getNumberOfTimeSteps,
           "Get how many advection steps were performed after the last apply() "
           "call.")
      .def("getTimeStepRatio", &Advect<T, D>::getTimeStepRatio,
           "Get the time step ratio used for advection.")
      .def("getCurrentTimeStep", &Advect<T, D>::getCurrentTimeStep,
           "Get the current time step.")
      .def("getCalculateNormalVectors",
           &Advect<T, D>::getCalculateNormalVectors,
           "Get whether normal vectors are computed during advection.")
      .def("setIntegrationScheme", &Advect<T, D>::setIntegrationScheme,
           "Set the integration scheme to use during advection.")
      .def("setDissipationAlpha", &Advect<T, D>::setDissipationAlpha,
           "Set the dissipation value to use for Lax Friedrichs integration.")
      .def("prepareLS", &Advect<T, D>::prepareLS, "Prepare the level-set.")
      // need scoped release since we are calling a python method from
      // parallelised C++ code here
      .def("apply", &Advect<T, D>::apply,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Perform advection.")
      .def("applyIntegration", &Advect<T, D>::applyIntegration,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Apply the integration scheme and calculate rates and maximum time "
           "step, but it do **not** move the surface.");

  // enums
  pybind11::enum_<IntegrationSchemeEnum>(module, "IntegrationSchemeEnum")
      .value("ENGQUIST_OSHER_1ST_ORDER",
             IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER)
      .value("ENGQUIST_OSHER_2ND_ORDER",
             IntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER)
      .value("LAX_FRIEDRICHS_1ST_ORDER",
             IntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER)
      .value("LAX_FRIEDRICHS_2ND_ORDER",
             IntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER",
             IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER)
      .value("LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
      .value("LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER",
             IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
      .value("LOCAL_LAX_FRIEDRICHS_2ND_ORDER",
             IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER)
      .value("STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER",
             IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);

  // BooleanOperation
  pybind11::class_<BooleanOperation<T, D>,
                   SmartPointer<BooleanOperation<T, D>>>(module,
                                                         "BooleanOperation")
      // constructors
      .def(pybind11::init(&SmartPointer<BooleanOperation<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<BooleanOperation<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<BooleanOperation<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init(
          [](SmartPointer<Domain<T, D>> &domain, BooleanOperationEnum op) {
            return SmartPointer<BooleanOperation<T, D>>::New(domain, op);
          }))
      .def(pybind11::init([](SmartPointer<Domain<T, D>> &domainA,
                             SmartPointer<Domain<T, D>> &domainB,
                             BooleanOperationEnum op) {
        return SmartPointer<BooleanOperation<T, D>>::New(domainA, domainB, op);
      }))
      // methods
      .def("setLevelset", &BooleanOperation<T, D>::setLevelSet,
           "Set levelset on which the boolean operation should be performed.")
      .def("setSecondLevelSet", &BooleanOperation<T, D>::setSecondLevelSet,
           "Set second levelset for boolean operation.")
      .def("setBooleanOperation", &BooleanOperation<T, D>::setBooleanOperation,
           "Set which type of boolean operation should be performed.")
      .def("apply", &BooleanOperation<T, D>::apply,
           "Perform the boolean operation.");
  // enums
  pybind11::enum_<BooleanOperationEnum>(module, "BooleanOperationEnum",
                                        pybind11::module_local())
      .value("INTERSECT", BooleanOperationEnum::INTERSECT)
      .value("UNION", BooleanOperationEnum::UNION)
      .value("RELATIVE_COMPLEMENT", BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .value("INVERT", BooleanOperationEnum::INVERT);

  pybind11::class_<CalculateCurvatures<T, D>,
                   SmartPointer<CalculateCurvatures<T, D>>>(
      module, "CalculateCurvatures")
      // constructors
      .def(pybind11::init(&SmartPointer<CalculateCurvatures<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<CalculateCurvatures<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init(
          [](SmartPointer<Domain<T, D>> &domain, CurvatureEnum type) {
            return SmartPointer<CalculateCurvatures<T, D>>::New(domain, type);
          }))
      // methods
      .def("setLevelSet", &CalculateCurvatures<T, D>::setLevelSet,
           "Set levelset for which to calculate the curvatures.")
      .def("setCurvatureType", &CalculateCurvatures<T, D>::setCurvatureType,
           "Set which method to use for calculation: Defaults to mean "
           "curvature.")
      .def("setMaxValue", &CalculateCurvatures<T, D>::setMaxValue,
           "Curvatures will be calculated for all LS values < maxValue.")
      .def("apply", &CalculateCurvatures<T, D>::apply,
           "Perform normal vector calculation.");

  // enums
  pybind11::enum_<CurvatureEnum>(module, "CurvatureEnum",
                                 pybind11::module_local())
      .value("MEAN_CURVATURE", CurvatureEnum::MEAN_CURVATURE)
      .value("GAUSSIAN_CURVATURE", CurvatureEnum::GAUSSIAN_CURVATURE)
      .value("MEAN_AND_GAUSSIAN_CURVATURE",
             CurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE);

  // CalculateNormalVectors
  pybind11::class_<CalculateNormalVectors<T, D>,
                   SmartPointer<CalculateNormalVectors<T, D>>>(
      module, "CalculateNormalVectors")
      // constructors
      .def(pybind11::init(&SmartPointer<CalculateNormalVectors<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<CalculateNormalVectors<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &CalculateNormalVectors<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &CalculateNormalVectors<T, D>::apply,
           "Perform normal vector calculation.");

  // CalculateVisibilities
  pybind11::class_<CalculateVisibilities<T, D>,
                   SmartPointer<CalculateVisibilities<T, D>>>(
      module, "CalculateVisibilities")
      .def(pybind11::init(
          &SmartPointer<CalculateVisibilities<T, D>>::New<
              SmartPointer<Domain<T, D>> &, const Vec3D<T> &, std::string>))
      .def("apply", &CalculateVisibilities<T, D>::apply);

  // Check
  pybind11::class_<Check<T, D>, SmartPointer<Check<T, D>>>(module, "Check")
      // constructors
      .def(pybind11::init(&SmartPointer<Check<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Check<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &Check<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &Check<T, D>::apply, "Perform check.");

#if VIENNALS_PYTHON_DIMENSION == 2
  // CompareArea
  pybind11::class_<CompareArea<T, D>, SmartPointer<CompareArea<T, D>>>(
      module, "CompareArea")
      // constructors
      .def(pybind11::init(&SmartPointer<CompareArea<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<CompareArea<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetTarget", &CompareArea<T, D>::setLevelSetTarget,
           "Sets the target level set.")
      .def("setLevelSetSample", &CompareArea<T, D>::setLevelSetSample,
           "Sets the sample level set.")
      .def("setDefaultIncrement", &CompareArea<T, D>::setDefaultIncrement,
           "Set default increment value")
      .def("setXRangeAndIncrement", &CompareArea<T, D>::setXRangeAndIncrement,
           "Sets the x-range and custom increment value")
      .def("setYRangeAndIncrement", &CompareArea<T, D>::setYRangeAndIncrement,
           "Sets the y-range and custom increment value")
      .def("setOutputMesh", &CompareArea<T, D>::setOutputMesh,
           "Set the output mesh where difference areas will be stored")
      .def("getAreaMismatch", &CompareArea<T, D>::getAreaMismatch,
           "Returns the computed area mismatch.")
      .def(
          "getCustomAreaMismatch", &CompareArea<T, D>::getCustomAreaMismatch,
          "Returns the computed area mismatch, with custom increments applied.")
      .def("getCellCount", &CompareArea<T, D>::getCellCount,
           "Returns the number of cells where the level sets differ.")
      .def("getCustomCellCount", &CompareArea<T, D>::getCustomCellCount,
           "Returns the number of cells where the level sets differ, with "
           "custom increments applied.")
      .def("apply", &CompareArea<T, D>::apply,
           "Computes the area difference between the two level sets.");

  // CompareNarrowBand
  pybind11::class_<CompareNarrowBand<T, D>,
                   SmartPointer<CompareNarrowBand<T, D>>>(module,
                                                          "CompareNarrowBand")
      // constructors
      .def(pybind11::init(&SmartPointer<CompareNarrowBand<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<CompareNarrowBand<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetTarget", &CompareNarrowBand<T, D>::setLevelSetTarget,
           "Sets the target level set.")
      .def("setLevelSetSample", &CompareNarrowBand<T, D>::setLevelSetSample,
           "Sets the sample level set.")
      .def("setXRange", &CompareNarrowBand<T, D>::setXRange,
           "Set the x-coordinate range to restrict the comparison area")
      .def("setYRange", &CompareNarrowBand<T, D>::setYRange,
           "Set the y-coordinate range to restrict the comparison area")
      .def("clearXRange", &CompareNarrowBand<T, D>::clearXRange,
           "Clear the x-range restriction")
      .def("clearYRange", &CompareNarrowBand<T, D>::clearYRange,
           "Clear the y-range restriction")
      .def("setOutputMesh", &CompareNarrowBand<T, D>::setOutputMesh,
           "Set the output mesh where difference values will be stored")
      .def("setOutputMeshSquaredDifferences",
           &CompareNarrowBand<T, D>::setOutputMeshSquaredDifferences,
           "Set whether to output squared differences (true) or absolute "
           "differences (false)")
      .def("apply", &CompareNarrowBand<T, D>::apply,
           "Apply the comparison and calculate the sum of squared differences.")
      .def("getSumSquaredDifferences",
           &CompareNarrowBand<T, D>::getSumSquaredDifferences,
           "Return the sum of squared differences calculated by apply().")
      .def("getSumDifferences", &CompareNarrowBand<T, D>::getSumDifferences,
           "Return the sum of absolute differences calculated by apply().")
      .def("getNumPoints", &CompareNarrowBand<T, D>::getNumPoints,
           "Return the number of points used in the comparison.")
      .def("getRMSE", &CompareNarrowBand<T, D>::getRMSE,
           "Calculate the root mean square error from previously computed "
           "values.");

  // CompareSparseField
  pybind11::class_<CompareSparseField<T, D>,
                   SmartPointer<CompareSparseField<T, D>>>(module,
                                                           "CompareSparseField")
      // constructors
      .def(pybind11::init(&SmartPointer<CompareSparseField<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<CompareSparseField<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetExpanded",
           &CompareSparseField<T, D>::setLevelSetExpanded,
           "Sets the expanded level set for comparison.")
      .def("setLevelSetIterated",
           &CompareSparseField<T, D>::setLevelSetIterated,
           "Sets the iterated level set to compare against the expanded one.")
      .def("setXRange", &CompareSparseField<T, D>::setXRange,
           "Set the x-coordinate range to restrict the comparison area")
      .def("setYRange", &CompareSparseField<T, D>::setYRange,
           "Set the y-coordinate range to restrict the comparison area")
      .def("clearXRange", &CompareSparseField<T, D>::clearXRange,
           "Clear the x-range restriction")
      .def("clearYRange", &CompareSparseField<T, D>::clearYRange,
           "Clear the y-range restriction")
      .def("setOutputMesh", &CompareSparseField<T, D>::setOutputMesh,
           "Set the output mesh where difference values will be stored")
      .def("setFillIteratedWithDistances",
           &CompareSparseField<T, D>::setFillIteratedWithDistances,
           "Set whether to fill the iterated level set with distance values")
      .def("setExpandedLevelSetWidth",
           &CompareSparseField<T, D>::setExpandedLevelSetWidth,
           "Set the expansion width for the expanded level set")
      .def("apply", &CompareSparseField<T, D>::apply,
           "Apply the comparison and calculate the sum of squared differences.")
      .def("getSumSquaredDifferences",
           &CompareSparseField<T, D>::getSumSquaredDifferences,
           "Return the sum of squared differences calculated by apply().")
      .def("getSumDifferences", &CompareSparseField<T, D>::getSumDifferences,
           "Return the sum of absolute differences calculated by apply().")
      .def("getNumPoints", &CompareSparseField<T, D>::getNumPoints,
           "Return the number of points used in the comparison.")
      .def("getNumSkippedPoints",
           &CompareSparseField<T, D>::getNumSkippedPoints,
           "Return the number of points skipped during comparison.")
      .def("getRMSE", &CompareSparseField<T, D>::getRMSE,
           "Calculate the root mean square error from previously computed "
           "values.");

  // CompareCriticalDimensions
  pybind11::class_<CompareCriticalDimensions<T, D>,
                   SmartPointer<CompareCriticalDimensions<T, D>>>(
      module, "CompareCriticalDimensions")
      // constructors
      .def(pybind11::init(&SmartPointer<CompareCriticalDimensions<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<CompareCriticalDimensions<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetReference",
           &CompareCriticalDimensions<T, D>::setLevelSetReference,
           "Sets the reference level set.")
      .def("setLevelSetCompare",
           &CompareCriticalDimensions<T, D>::setLevelSetCompare,
           "Sets the comparison level set.")
      .def("addXRange", &CompareCriticalDimensions<T, D>::addXRange,
           pybind11::arg("minX"), pybind11::arg("maxX"),
           pybind11::arg("findMaximum") = true,
           "Add an X range to find maximum or minimum Y position.")
      .def("addYRange", &CompareCriticalDimensions<T, D>::addYRange,
           pybind11::arg("minY"), pybind11::arg("maxY"),
           pybind11::arg("findMaximum") = true,
           "Add a Y range to find maximum or minimum X position.")
      .def("clearRanges", &CompareCriticalDimensions<T, D>::clearRanges,
           "Clear all range specifications.")
      .def("setOutputMesh", &CompareCriticalDimensions<T, D>::setOutputMesh,
           "Set the output mesh where critical dimension locations will be "
           "stored.")
      .def("apply", &CompareCriticalDimensions<T, D>::apply,
           "Apply the comparison.")
      .def("getNumCriticalDimensions",
           &CompareCriticalDimensions<T, D>::getNumCriticalDimensions,
           "Get the number of critical dimensions compared.")
      .def("getCriticalDimensionResult",
           [](CompareCriticalDimensions<T, D> &self, size_t index) {
             T posRef, posCmp, diff;
             bool valid = self.getCriticalDimensionResult(index, posRef, posCmp,
                                                          diff);
             if (valid) {
               return pybind11::make_tuple(true, posRef, posCmp, diff);
             } else {
               return pybind11::make_tuple(false, 0.0, 0.0, 0.0);
             }
           },
           pybind11::arg("index"),
           "Get a specific critical dimension result. Returns (valid, "
           "positionReference, positionCompare, difference).")
      .def("getMeanDifference",
           &CompareCriticalDimensions<T, D>::getMeanDifference,
           "Get mean absolute difference across all valid critical dimensions.")
      .def("getMaxDifference", &CompareCriticalDimensions<T, D>::getMaxDifference,
           "Get maximum difference across all valid critical dimensions.")
      .def("getRMSE", &CompareCriticalDimensions<T, D>::getRMSE,
           "Get RMSE across all valid critical dimensions.")
      .def("getAllDifferences",
           &CompareCriticalDimensions<T, D>::getAllDifferences,
           "Get all valid differences as a list.");

  // CompareChamfer
  pybind11::class_<CompareChamfer<T, D>, SmartPointer<CompareChamfer<T, D>>>(
      module, "CompareChamfer")
      // constructors
      .def(pybind11::init(&SmartPointer<CompareChamfer<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<CompareChamfer<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                   SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetTarget", &CompareChamfer<T, D>::setLevelSetTarget,
           "Set the target level set.")
      .def("setLevelSetSample", &CompareChamfer<T, D>::setLevelSetSample,
           "Set the sample level set.")
      .def("setOutputMeshTarget", &CompareChamfer<T, D>::setOutputMeshTarget,
           "Set output mesh for target surface points with distance data.")
      .def("setOutputMeshSample", &CompareChamfer<T, D>::setOutputMeshSample,
           "Set output mesh for sample surface points with distance data.")
      .def("apply", &CompareChamfer<T, D>::apply,
           "Apply the Chamfer distance calculation.")
      .def("getForwardDistance", &CompareChamfer<T, D>::getForwardDistance,
           "Get the forward distance (average distance from target to sample).")
      .def("getBackwardDistance", &CompareChamfer<T, D>::getBackwardDistance,
           "Get the backward distance (average distance from sample to target).")
      .def("getChamferDistance", &CompareChamfer<T, D>::getChamferDistance,
           "Get the Chamfer distance (average of forward and backward).")
      .def("getRMSChamferDistance",
           &CompareChamfer<T, D>::getRMSChamferDistance,
           "Get the RMS Chamfer distance.")
      .def("getMaxDistance", &CompareChamfer<T, D>::getMaxDistance,
           "Get the maximum nearest-neighbor distance.")
      .def("getNumTargetPoints", &CompareChamfer<T, D>::getNumTargetPoints,
           "Get the number of target surface points.")
      .def("getNumSamplePoints", &CompareChamfer<T, D>::getNumSamplePoints,
           "Get the number of sample surface points.");
#endif

  // ConvexHull
  pybind11::class_<ConvexHull<T, D>, SmartPointer<ConvexHull<T, D>>>(
      module, "ConvexHull")
      // constructors
      .def(pybind11::init(&SmartPointer<ConvexHull<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<ConvexHull<T, D>>::New<
              SmartPointer<Mesh<T>> &, SmartPointer<PointCloud<T, D>> &>))
      // methods
      .def("setMesh", &ConvexHull<T, D>::setMesh,
           "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &ConvexHull<T, D>::setPointCloud,
           "Set point cloud used to generate mesh.")
      .def("apply", &ConvexHull<T, D>::apply, "Generate Hull.");

  // DetectFeatures
  pybind11::class_<DetectFeatures<T, D>, SmartPointer<DetectFeatures<T, D>>>(
      module, "DetectFeatures")
      // constructors
      .def(pybind11::init(&SmartPointer<DetectFeatures<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<DetectFeatures<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<DetectFeatures<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                   T>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(pybind11::init([](SmartPointer<Domain<T, D>> &domain, T maxValue,
                             FeatureDetectionEnum type) {
        return SmartPointer<DetectFeatures<T, D>>::New(domain, maxValue, type);
      }))
      .def("setDetectionThreshold",
           &DetectFeatures<T, D>::setDetectionThreshold,
           "Set the curvature value above which a point is considered a "
           "feature.")
      .def("setDetectionMethod", &DetectFeatures<T, D>::setDetectionMethod,
           "Set which method to use to detect features. Defaults to Curvature.")
      .def("apply", &DetectFeatures<T, D>::apply, "Detect features.");

  // enums
  pybind11::enum_<FeatureDetectionEnum>(module, "FeatureDetectionEnum",
                                        pybind11::module_local())
      .value("CURVATURE", FeatureDetectionEnum::CURVATURE)
      .value("NORMALS_ANGLE", FeatureDetectionEnum::NORMALS_ANGLE);

  // Domain
  pybind11::class_<Domain<T, D>, SmartPointer<Domain<T, D>>>(module, "Domain")
      // constructors
      .def(pybind11::init(&SmartPointer<Domain<T, D>>::New<>))
      .def(pybind11::init(
               &SmartPointer<Domain<T, D>>::New<viennahrle::CoordType>),
           pybind11::arg("gridDelta") = 1.0)
      //  .def(pybind11::init(
      //      &SmartPointer<Domain<T, D>>::New<viennahrle::CoordType *,
      //                                       BoundaryConditionEnum *>))
      .def(pybind11::init([](std::array<viennahrle::CoordType, 2 * D> bounds,
                             std::array<BoundaryConditionEnum, D> bcs,
                             viennahrle::CoordType gridDelta) {
             return SmartPointer<Domain<T, D>>::New(bounds.data(), bcs.data(),
                                                    gridDelta);
           }),
           pybind11::arg("bounds"), pybind11::arg("boundaryConditions"),
           pybind11::arg("gridDelta") = 1.0)
      .def(pybind11::init(&SmartPointer<Domain<T, D>>::New<
                          std::vector<viennahrle::CoordType>,
                          std::vector<unsigned>, viennahrle::CoordType>),
           pybind11::arg("bounds"), pybind11::arg("boundaryConditions"),
           pybind11::arg("gridDelta") = 1.0)
      //  .def(pybind11::init(
      //      &SmartPointer<Domain<T, D>>::New<Domain<T,
      //      D>::PointValueVectorType,
      //                                       viennahrle::CoordType *,
      //                                       BoundaryConditionEnum *>))
      //  .def(pybind11::init(&SmartPointer<Domain<T, D>>::New<
      //                      Domain<T, D>::PointValueVectorType,
      //                      viennahrle::CoordType
      //                      *, BoundaryConditionEnum *,
      //                      viennahrle::CoordType>))
      .def(pybind11::init(
          &SmartPointer<Domain<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Domain<T, D>>::New<viennahrle::Grid<D> &>))
      // methods
      .def("deepCopy", &Domain<T, D>::deepCopy,
           "Copy lsDomain in this lsDomain.")
      .def("getNumberOfSegments", &Domain<T, D>::getNumberOfSegments,
           "Get the number of segments, the level set structure is divided "
           "into.")
      .def("getNumberOfPoints", &Domain<T, D>::getNumberOfPoints,
           "Get the number of defined level set values.")
      .def("getLevelSetWidth", &Domain<T, D>::getLevelSetWidth,
           "Get the number of layers of level set points around the explicit "
           "surface.")
      .def("setLevelSetWidth", &Domain<T, D>::setLevelSetWidth,
           "Set the number of layers of level set points which should be "
           "stored around the explicit surface.")
      .def("clearMetaData", &Domain<T, D>::clearMetaData,
           "Clear all metadata stored in the level set.")
      // allow filehandle to be passed and default to python standard output
      .def(
          "print",
          [](Domain<T, D> &d, pybind11::object fileHandle) {
            if (!(pybind11::hasattr(fileHandle, "write") &&
                  pybind11::hasattr(fileHandle, "flush"))) {
              throw pybind11::type_error(
                  "MyClass::read_from_file_like_object(file): incompatible "
                  "function argument:  `file` must be a file-like object, but "
                  "`" +
                  (std::string)(pybind11::repr(fileHandle)) + "` provided");
            }
            pybind11::detail::pythonbuf buf(fileHandle);
            std::ostream stream(&buf);
            d.print(stream);
          },
          pybind11::arg("stream") =
              pybind11::module::import("sys").attr("stdout"));

  // enums
  pybind11::enum_<BoundaryConditionEnum>(module, "BoundaryConditionEnum")
      .value("REFLECTIVE_BOUNDARY", BoundaryConditionEnum::REFLECTIVE_BOUNDARY)
      .value("INFINITE_BOUNDARY", BoundaryConditionEnum::INFINITE_BOUNDARY)
      .value("PERIODIC_BOUNDARY", BoundaryConditionEnum::PERIODIC_BOUNDARY)
      .value("POS_INFINITE_BOUNDARY",
             BoundaryConditionEnum::POS_INFINITE_BOUNDARY)
      .value("NEG_INFINITE_BOUNDARY",
             BoundaryConditionEnum::NEG_INFINITE_BOUNDARY);

  // GeometricAdvect
  pybind11::class_<GeometricAdvect<T, D>, SmartPointer<GeometricAdvect<T, D>>>(
      module, "GeometricAdvect")
      // constructors
      .def(pybind11::init(&SmartPointer<GeometricAdvect<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<GeometricAdvect<T, D>>::New<
              SmartPointer<Domain<T, D>> &,
              SmartPointer<
                  GeometricAdvectDistribution<viennahrle::CoordType, D>> &>))
      // methods
      .def("setLevelSet", &GeometricAdvect<T, D>::setLevelSet,
           "Set levelset to advect.")
      .def(
          "setAdvectionDistribution",
          &GeometricAdvect<T, D>::setAdvectionDistribution<
              PylsGeometricAdvectDistribution>,
          "Set advection distribution to use as kernel for the fast advection.")
      .def("apply", &GeometricAdvect<T, D>::apply,
           pybind11::call_guard<pybind11::gil_scoped_release>(),
           "Perform advection.");

  // GeometricAdvectDistributions
  pybind11::class_<GeometricAdvectDistribution<T, D>,
                   SmartPointer<GeometricAdvectDistribution<T, D>>,
                   PylsGeometricAdvectDistribution>(
      module, "GeometricAdvectDistribution")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("isInside", &GeometricAdvectDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance",
           &GeometricAdvectDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &GeometricAdvectDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<SphereDistribution<T, D>,
                   SmartPointer<SphereDistribution<T, D>>,
                   GeometricAdvectDistribution<T, D>>(module,
                                                      "SphereDistribution")
      // constructors
      .def(pybind11::init(&SmartPointer<SphereDistribution<T, D>>::New<T, T>))
      // methods
      .def("isInside", &SphereDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &SphereDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &SphereDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  pybind11::class_<BoxDistribution<T, D>, SmartPointer<BoxDistribution<T, D>>,
                   GeometricAdvectDistribution<T, D>>(module, "BoxDistribution")
      // constructors
      .def(pybind11::init(
          &SmartPointer<BoxDistribution<T, D>>::New<const std::array<T, 3>, T>))
      // methods
      .def("isInside", &BoxDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &BoxDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &BoxDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  // lsExpand
  pybind11::class_<Expand<T, D>, SmartPointer<Expand<T, D>>>(module, "Expand")
      // constructors
      .def(pybind11::init(&SmartPointer<Expand<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Expand<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Expand<T, D>>::New<SmartPointer<Domain<T, D>> &, int>))
      // methods
      .def("setLevelSet", &Expand<T, D>::setLevelSet, "Set levelset to expand.")
      .def("setWidth", &Expand<T, D>::setWidth, "Set the width to expand to.")
      .def("apply", &Expand<T, D>::apply, "Perform expansion.");

  // lsExtrude
  // Does not work in current implementation, because one can not import both 2D
  // and 3D ViennaLS libraries in Python in the same file
  //   pybind11::class_<lsExtrude<T>, SmartPointer<lsExtrude<T>>>(module,
  //                                                                "lsExtrude")
  //       // constructors
  //       .def(pybind11::init(&SmartPointer<lsExtrude<T>>::New<>))
  //       .def(
  //           pybind11::init(&SmartPointer<lsExtrude<T>>::New<
  //                           SmartPointer< Domain<T, 2>> &,
  //                           SmartPointer< Domain<T, 3>> &,
  //                          std::array<T, 2>, const int,
  //                          std::array<BoundaryConditionEnum<3>, 3>>))
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
  //            pybind11::overload_cast<std::array<BoundaryConditionEnum<3>,
  //            3>>(
  //                &lsExtrude<T>::setBoundaryConditions),
  //            "Set the boundary conditions in the 3D extruded domain.")
  //       .def("setBoundaryConditions",
  //            pybind11::overload_cast<BoundaryConditionEnum<3> *>(
  //                &lsExtrude<T>::setBoundaryConditions),
  //            "Set the boundary conditions in the 3D extruded domain.")
  //       .def("apply", &lsExtrude<T>::apply, "Perform extrusion.");

  // lsFileFormats
  pybind11::enum_<FileFormatEnum>(module, "FileFormatEnum",
                                  pybind11::module_local())
      .value("VTK_LEGACY", FileFormatEnum::VTK_LEGACY)
      .value("VTP", FileFormatEnum::VTP)
      .value("VTU", FileFormatEnum::VTU);

  // FromSurfaceMesh
  pybind11::class_<FromSurfaceMesh<T, D>, SmartPointer<FromSurfaceMesh<T, D>>>(
      module, "FromSurfaceMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<FromSurfaceMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<FromSurfaceMesh<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(
          &SmartPointer<FromSurfaceMesh<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &, bool>))
      // methods
      .def("setLevelSet", &FromSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to read into.")
      .def("setMesh", &FromSurfaceMesh<T, D>::setMesh,
           "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles",
           static_cast<void (FromSurfaceMesh<T, D>::*)(bool)>(
               &FromSurfaceMesh<T, D>::setRemoveBoundaryTriangles),
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("setRemoveBoundaryTriangles",
           static_cast<void (FromSurfaceMesh<T, D>::*)(std::array<bool, 3>)>(
               &FromSurfaceMesh<T, D>::setRemoveBoundaryTriangles),
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("apply", &FromSurfaceMesh<T, D>::apply,
           "Construct a levelset from a surface mesh.");

  // FromVolumeMesh
  pybind11::class_<FromVolumeMesh<T, D>, SmartPointer<FromVolumeMesh<T, D>>>(
      module, "FromVolumeMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<FromVolumeMesh<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<FromVolumeMesh<T, D>>::New<
                          typename FromVolumeMesh<T, D>::GridType &,
                          SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(&SmartPointer<FromVolumeMesh<T, D>>::New<
                          typename FromVolumeMesh<T, D>::GridType &,
                          SmartPointer<Mesh<T>> &, bool>))
      // methods
      .def("setGrid", &FromVolumeMesh<T, D>::setGrid,
           "Set the grid used to read in the level sets.")
      .def("setMesh", &FromVolumeMesh<T, D>::setMesh,
           "Set the mesh to read from.")
      .def("setRemoveBoundaryTriangles",
           &FromVolumeMesh<T, D>::setRemoveBoundaryTriangles,
           "Set whether to include mesh elements outside of the simulation "
           "domain.")
      .def("apply", &FromVolumeMesh<T, D>::apply,
           "Construct a levelset from a volume mesh.");

  // FromMesh
  pybind11::class_<FromMesh<T, D>, SmartPointer<FromMesh<T, D>>>(module,
                                                                 "FromMesh")
      .def(pybind11::init(&SmartPointer<FromMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<FromMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                             SmartPointer<Mesh<T>> &>))
      .def("setMesh", &FromMesh<T, D>::setMesh, "Set the mesh to read from.")
      .def("setSortPointList", &FromMesh<T, D>::setSortPointList)
      .def("apply", &FromMesh<T, D>::apply);

  // lsGeometries
  // Sphere
  pybind11::class_<Sphere<T, D>, SmartPointer<Sphere<T, D>>>(module, "Sphere")
      // constructors
      .def(pybind11::init(&SmartPointer<Sphere<T, D>>::New<
                          const std::vector<T> & /*origin*/, T /*radius*/>),
           pybind11::arg("origin"), pybind11::arg("radius"));
  // Plane
  pybind11::class_<Plane<T, D>, SmartPointer<Plane<T, D>>>(module, "Plane")
      // constructors
      .def(pybind11::init(&SmartPointer<Plane<T, D>>::New<
                          const std::vector<T> & /*origin*/,
                          const std::vector<T> & /*normal*/>),
           pybind11::arg("origin"), pybind11::arg("normal"));
  // Box
  pybind11::class_<Box<T, D>, SmartPointer<Box<T, D>>>(module, "Box")
      // constructors
      .def(pybind11::init(&SmartPointer<Box<T, D>>::New<
                          const std::vector<T> & /*minPoint*/,
                          const std::vector<T> & /*maxPoint*/>),
           pybind11::arg("minPoint"), pybind11::arg("maxPoint"));
  // Cylinder
  pybind11::class_<Cylinder<T, D>, SmartPointer<Cylinder<T, D>>>(module,
                                                                 "Cylinder")
      // constructors
      .def(pybind11::init(
               &SmartPointer<Cylinder<T, D>>::New<
                   const std::vector<T> & /*origin*/,
                   const std::vector<T> & /*axisDirection*/, const T /*height*/,
                   const T /*radius*/, const T /*topRadius*/>),
           pybind11::arg("origin"), pybind11::arg("axisDirection"),
           pybind11::arg("height"), pybind11::arg("radius"),
           pybind11::arg("topRadius") = 0.);
  // PointCloud
  pybind11::class_<PointCloud<T, D>, SmartPointer<PointCloud<T, D>>>(
      module, "PointCloud")
      // constructors
      .def(pybind11::init(&SmartPointer<PointCloud<T, D>>::New<
                          const std::vector<VectorType<T, D>> &>))
      // methods
      .def("insertNextPoint",
           (void(PointCloud<T, D>::*)(const VectorType<T, D> &)) &
               PointCloud<T, D>::insertNextPoint);

  // MakeGeometry
  pybind11::class_<MakeGeometry<T, D>, SmartPointer<MakeGeometry<T, D>>>(
      module, "MakeGeometry")
      // constructors
      .def(pybind11::init(&SmartPointer<MakeGeometry<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                 SmartPointer<Sphere<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                 SmartPointer<Plane<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                 SmartPointer<Box<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Cylinder<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MakeGeometry<T, D>>::New<
              SmartPointer<Domain<T, D>> &, SmartPointer<PointCloud<T, D>> &>))
      // methods
      .def("setLevelSet", &MakeGeometry<T, D>::setLevelSet,
           "Set the levelset in which to create the geometry.")
      .def("setGeometry",
           (void(MakeGeometry<T, D>::*)(SmartPointer<Sphere<T, D>>)) &
               MakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(MakeGeometry<T, D>::*)(SmartPointer<Plane<T, D>>)) &
               MakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(MakeGeometry<T, D>::*)(SmartPointer<Box<T, D>>)) &
               MakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(MakeGeometry<T, D>::*)(SmartPointer<Cylinder<T, D>>)) &
               MakeGeometry<T, D>::setGeometry)
      .def("setGeometry",
           (void(MakeGeometry<T, D>::*)(SmartPointer<PointCloud<T, D>>)) &
               MakeGeometry<T, D>::setGeometry)
      .def("setIgnoreBoundaryConditions",
           (void(MakeGeometry<T, D>::*)(bool)) &
               MakeGeometry<T, D>::setIgnoreBoundaryConditions)
      .def("setIgnoreBoundaryConditions",
           (void(MakeGeometry<T, D>::*)(std::array<bool, 3>)) &
               MakeGeometry<T, D>::setIgnoreBoundaryConditions)
      .def("apply", &MakeGeometry<T, D>::apply, "Generate the geometry.");

  // MarkVoidPoints
  pybind11::class_<MarkVoidPoints<T, D>, SmartPointer<MarkVoidPoints<T, D>>>(
      module, "MarkVoidPoints")
      // constructors
      .def(pybind11::init(&SmartPointer<MarkVoidPoints<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<MarkVoidPoints<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<MarkVoidPoints<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                   bool &>))
      // methods
      .def("setLevelSet", &MarkVoidPoints<T, D>::setLevelSet,
           "Set the levelset to mark void points in.")
      .def("setReverseVoidDetection",
           &MarkVoidPoints<T, D>::setReverseVoidDetection,
           "Reverse the logic of detecting the top surface.")
      .def("setDetectLargestSurface",
           &MarkVoidPoints<T, D>::setDetectLargestSurface,
           "Set that the top surface should be the one with the most connected "
           "LS points.")
      .def("setVoidTopSurface", &MarkVoidPoints<T, D>::setVoidTopSurface,
           "Set the logic by which to choose the surface which is non-void. "
           "All other connected surfaces will then be marked as void points.")
      .def("setSaveComponentsId", &MarkVoidPoints<T, D>::setSaveComponentIds,
           "Save the connectivity information of all LS points in the "
           "pointData of the level set.")
      .def("apply", &MarkVoidPoints<T, D>::apply, "Mark void points.");

  // VoidTopSurfaceEnum
  pybind11::enum_<VoidTopSurfaceEnum>(module, "VoidTopSurfaceEnum",
                                      pybind11::module_local())
      .value("LEX_LOWEST", VoidTopSurfaceEnum::LEX_LOWEST)
      .value("LEX_HIGHEST", VoidTopSurfaceEnum::LEX_HIGHEST)
      .value("LARGEST", VoidTopSurfaceEnum::LARGEST)
      .value("SMALLEST", VoidTopSurfaceEnum::SMALLEST);

  // PointData
  pybind11::class_<PointData<T>, SmartPointer<PointData<T>>>(module,
                                                             "PointData")
      // constructors
      .def(pybind11::init(&SmartPointer<PointData<T>>::New<>))
      // methods
      .def("insertNextScalarData",
           (void(PointData<T>::*)(const PointData<T>::ScalarDataType &,
                                  const std::string &)) &
               PointData<T>::insertNextScalarData,
           pybind11::arg("scalars"), pybind11::arg("label") = "Scalars")
      .def("insertNextVectorData",
           (void(PointData<T>::*)(const PointData<T>::VectorDataType &,
                                  const std::string &)) &
               PointData<T>::insertNextVectorData,
           pybind11::arg("vectors"), pybind11::arg("label") = "Vectors")
      .def("getScalarDataSize", &PointData<T>::getScalarDataSize)
      .def("getVectorDataSize", &PointData<T>::getVectorDataSize)
      .def("getScalarData",
           (PointData<T>::ScalarDataType * (PointData<T>::*)(int)) &
               PointData<T>::getScalarData)
      .def("getScalarData", (PointData<T>::ScalarDataType *
                             (PointData<T>::*)(const std::string &, bool)) &
                                PointData<T>::getScalarData)
      .def("getScalarDataLabel", &PointData<T>::getScalarDataLabel)
      .def("getVectorData",
           (PointData<T>::VectorDataType * (PointData<T>::*)(int)) &
               PointData<T>::getVectorData)
      .def("getVectorData", (PointData<T>::VectorDataType *
                             (PointData<T>::*)(const std::string &, bool)) &
                                PointData<T>::getVectorData)
      .def("getVectorDataLabel", &PointData<T>::getVectorDataLabel);

  // MaterialMap
  pybind11::class_<MaterialMap, SmartPointer<MaterialMap>>(module,
                                                           "MaterialMap")
      // constructors
      .def(pybind11::init(&SmartPointer<MaterialMap>::New<>))
      // methods
      .def("insertNextMaterial", &MaterialMap::insertNextMaterial,
           "Insert a new material into the map.")
      .def("setMaterialId", &MaterialMap::setMaterialId)
      .def("getNumberOfLayers", &MaterialMap::getNumberOfLayers,
           "Get the number of level-sets in the material map.")
      .def("getNumberOfMaterials", &MaterialMap::getNumberOfMaterials)
      .def("getMaterialId", &MaterialMap::getMaterialId);

  // Mesh<T>
  pybind11::class_<Mesh<T>, SmartPointer<Mesh<T>>>(module, "Mesh")
      // constructors
      .def(pybind11::init(&SmartPointer<Mesh<T>>::New<>))
      // methods
      .def("getNodes",
           (std::vector<std::array<T, 3>> & (Mesh<T>::*)()) & Mesh<T>::getNodes,
           "Get all nodes of the mesh as a list.")
      .def("getVerticies",
           (std::vector<std::array<unsigned, 1>> & (Mesh<T>::*)()) &
               Mesh<T>::getElements<1>,
           "Get a list of verticies of the mesh.")
      .def("getLines",
           (std::vector<std::array<unsigned, 2>> & (Mesh<T>::*)()) &
               Mesh<T>::getElements<2>,
           "Get a list of lines of the mesh.")
      .def("getTriangles",
           (std::vector<std::array<unsigned, 3>> & (Mesh<T>::*)()) &
               Mesh<T>::getElements<3>,
           "Get a list of triangles of the mesh.")
      .def("getTetras",
           (std::vector<std::array<unsigned, 4>> & (Mesh<T>::*)()) &
               Mesh<T>::getElements<4>,
           "Get a list of tetrahedrons of the mesh.")
      .def("getHexas",
           (std::vector<std::array<unsigned, 8>> & (Mesh<T>::*)()) &
               Mesh<T>::getElements<8>,
           "Get a list of hexahedrons of the mesh.")
      .def("getPointData",
           (PointData<T> & (Mesh<T>::*)()) & Mesh<T>::getPointData,
           "Return a reference to the point data of the mesh.")
      .def("getCellData",
           (PointData<T> & (Mesh<T>::*)()) & Mesh<T>::getCellData,
           "Return a reference to the cell data of the mesh.")
      .def("insertNextNode", &Mesh<T>::insertNextNode,
           "Insert a node in the mesh.")
      .def("insertNextVertex", &Mesh<T>::insertNextVertex,
           "Insert a vertex in the mesh.")
      .def("insertNextLine", &Mesh<T>::insertNextLine,
           "Insert a line in the mesh.")
      .def("insertNextTriangle", &Mesh<T>::insertNextTriangle,
           "Insert a triangle in the mesh.")
      .def("insertNextTetra", &Mesh<T>::insertNextTetra,
           "Insert a tetrahedron in the mesh.")
      .def("insertNextHexa", &Mesh<T>::insertNextHexa,
           "Insert a hexahedron in the mesh.")
      .def("removeDuplicateNodes", &Mesh<T>::removeDuplicateNodes,
           "Remove nodes which occur twice in the mesh, and replace their IDs "
           "in the mesh elements.")
      .def("append", &Mesh<T>::append, "Append another mesh to this mesh.")
      .def("print", &Mesh<T>::print, "Print basic information about the mesh.")
      .def("clear", &Mesh<T>::clear, "Clear all data in the mesh.");

  // TransformEnum
  pybind11::enum_<TransformEnum>(module, "TransformEnum",
                                 pybind11::module_local())
      .value("TRANSLATION", TransformEnum::TRANSLATION)
      .value("ROTATION", TransformEnum::ROTATION)
      .value("SCALE", TransformEnum::SCALE);

  // TransformMesh
  pybind11::class_<TransformMesh<T>, SmartPointer<TransformMesh<T>>>(
      module, "TransformMesh")
      // constructors
      .def(pybind11::init([](SmartPointer<Mesh<T>> &mesh, TransformEnum op,
                             Vec3D<T> vec, double angle) {
             return SmartPointer<TransformMesh<T>>::New(mesh, op, vec, angle);
           }),
           pybind11::arg("mesh"),
           pybind11::arg("transform") = TransformEnum::TRANSLATION,
           pybind11::arg("transformVector") = Vec3D<T>{0., 0., 0.},
           pybind11::arg("angle") = 0.)
      // methods
      .def("apply", &TransformMesh<T>::apply, "Apply the transformation.");

  // Prune
  pybind11::class_<Prune<T, D>, SmartPointer<Prune<T, D>>>(module, "Prune")
      // constructors
      .def(pybind11::init(&SmartPointer<Prune<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Prune<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &Prune<T, D>::setLevelSet, "Set levelset to prune.")
      .def("apply", &Prune<T, D>::apply, "Perform pruning operation.");

  // Reader
  pybind11::class_<Reader<T, D>, SmartPointer<Reader<T, D>>>(module, "Reader")
      // constructors
      .def(pybind11::init(&SmartPointer<Reader<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Reader<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Reader<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           std::string>))
      // methods
      .def("setLevelSet", &Reader<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &Reader<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &Reader<T, D>::apply, "Write to file.");

  // Reduce
  pybind11::class_<Reduce<T, D>, SmartPointer<Reduce<T, D>>>(module, "Reduce")
      // constructors
      .def(pybind11::init(&SmartPointer<Reduce<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Reduce<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Reduce<T, D>>::New<SmartPointer<Domain<T, D>> &, int>))
      .def(pybind11::init(
          &SmartPointer<Reduce<T, D>>::New<SmartPointer<Domain<T, D>> &, int,
                                           bool>))
      // methods
      .def("setLevelSet", &Reduce<T, D>::setLevelSet, "Set levelset to reduce.")
      .def("setWidth", &Reduce<T, D>::setWidth, "Set the width to reduce to.")
      .def("setNoNewSegment", &Reduce<T, D>::setNoNewSegment,
           "Set whether the levelset should be segmented anew (balanced across "
           "cores) after reduction.")
      .def("apply", &Reduce<T, D>::apply, "Perform reduction.");

  // RemoveStrayPoints
  pybind11::class_<RemoveStrayPoints<T, D>,
                   SmartPointer<RemoveStrayPoints<T, D>>>(module,
                                                          "RemoveStrayPoints")
      // constructors
      .def(pybind11::init(&SmartPointer<RemoveStrayPoints<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<RemoveStrayPoints<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &RemoveStrayPoints<T, D>::setLevelSet,
           "Set levelset for stray point removal.")
      .def("setVoidTopSurface", &RemoveStrayPoints<T, D>::setVoidTopSurface,
           "Set the logic by which to choose the surface which should be kept. "
           "All other LS values will be marked as stray points and removed.")
      .def("apply", &RemoveStrayPoints<T, D>::apply, "Remove stray points.");

  // Slice
  pybind11::class_<Slice<T>, SmartPointer<Slice<T>>>(module, "Slice")
      // constructors
      .def(pybind11::init(&SmartPointer<Slice<T>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Slice<T>>::New<SmartPointer<Domain<T, 3>> &,
                                       SmartPointer<Domain<T, 2>> &, int, T>))
      .def(pybind11::init(
          &SmartPointer<Slice<T>>::New<SmartPointer<Domain<T, 3>> &, int, T>))
      // methods
      .def("setSourceLevelSet", &Slice<T>::setSourceLevelSet,
           "Set the 3D source level set from which to extract the slice.")
      .def("setSliceLevelSet", &Slice<T>::setSliceLevelSet,
           "Set the 2D level set where the extracted slice will be stored.")
      .def("setSliceDimension", &Slice<T>::setSliceDimension,
           "Set the dimension along which to slice (0=x, 1=y, 2=z).")
      .def("setSlicePosition", &Slice<T>::setSlicePosition,
           "Set the position along the slice dimension where to extract the "
           "slice.")
      .def("setWritePath", &Slice<T>::setWritePath,
           "Set the path where the slice should be written to.")
      .def("getSliceLevelSet", &Slice<T>::getSliceLevelSet,
           "Get the 2D slice level set after extraction.")
      .def("setReflectX", &Slice<T>::setReflectX,
           "Set whether to reflect all x-coordinates in the resulting slice.")
      .def("apply", &Slice<T>::apply,
           "Extract the 2D slice from the 3D domain.");

  // ToDiskMesh
  pybind11::class_<ToDiskMesh<T, D>, SmartPointer<ToDiskMesh<T, D>>>(
      module, "ToDiskMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<ToDiskMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<ToDiskMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                               SmartPointer<Mesh<T>> &>))
      // methods
      .def("setLevelSet", &ToDiskMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &ToDiskMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &ToDiskMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToMesh
  pybind11::class_<ToMesh<T, D>, SmartPointer<ToMesh<T, D>>>(module, "ToMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<ToMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<ToMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(
          &SmartPointer<ToMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           SmartPointer<Mesh<T>> &, bool>))
      .def(pybind11::init(
          &SmartPointer<ToMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           SmartPointer<Mesh<T>> &, bool,
                                           bool>))
      // methods
      .def("setLevelSet", &ToMesh<T, D>::setLevelSet, "Set levelset to mesh.")
      .def("setMesh", &ToMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("setOnlyDefined", &ToMesh<T, D>::setOnlyDefined,
           "Set whether only defined points should be output to the mesh.")
      .def("setOnlyActive", &ToMesh<T, D>::setOnlyActive,
           "Set whether only level set points <0.5 should be output.")
      .def("apply", &ToMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToSurfaceMesh
  pybind11::class_<ToSurfaceMesh<T, D>, SmartPointer<ToSurfaceMesh<T, D>>>(
      module, "ToSurfaceMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<ToSurfaceMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<ToSurfaceMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                  SmartPointer<Mesh<T>> &>))
      // methods
      .def("setLevelSet", &ToSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &ToSurfaceMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("apply", &ToSurfaceMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToVoxelMesh
  pybind11::class_<ToVoxelMesh<T, D>, SmartPointer<ToVoxelMesh<T, D>>>(
      module, "ToVoxelMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<ToVoxelMesh<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<ToVoxelMesh<T, D>>::New<SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(
          &SmartPointer<ToVoxelMesh<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                                SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(&SmartPointer<ToVoxelMesh<T, D>>::New<
                          std::vector<SmartPointer<Domain<T, D>>> &,
                          SmartPointer<Mesh<T>> &>))
      // methods
      .def("insertNextLevelSet", &ToVoxelMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the mesh.")
      .def("setMesh", &ToVoxelMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &ToVoxelMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // VelocityField
  pybind11::class_<VelocityField<T>, SmartPointer<VelocityField<T>>,
                   PylsVelocityField>(module, "VelocityField")
      // constructors
      .def(pybind11::init<>())
      // methods
      .def("getScalarVelocity", &VelocityField<T>::getScalarVelocity,
           "Return the scalar velocity for a point of material at coordinate "
           "with normal vector normal.")
      .def("getVectorVelocity", &VelocityField<T>::getVectorVelocity,
           "Return the vector velocity for a point of material at coordinate "
           "with normal vector normal.")
      .def("getDissipationAlpha", &VelocityField<T>::getDissipationAlpha,
           "Return the analytical dissipation alpha value if the "
           "lsLocalLaxFriedrichsAnalytical scheme is used for advection.");

  // VTKReader
  pybind11::class_<VTKReader<T>, SmartPointer<VTKReader<T>>>(module,
                                                             "VTKReader")
      // constructors
      .def(pybind11::init(&SmartPointer<VTKReader<T>>::New<>))
      .def(pybind11::init(
          &SmartPointer<VTKReader<T>>::New<SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(
          &SmartPointer<VTKReader<T>>::New<SmartPointer<Mesh<T>> &,
                                           std::string>))
      .def(pybind11::init([](SmartPointer<Mesh<T>> &mesh, FileFormatEnum format,
                             std::string s) {
        return SmartPointer<VTKReader<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &VTKReader<T>::setMesh, "Set the mesh to read into.")
      .def("setFileFormat", &VTKReader<T>::setFileFormat,
           "Set the file format of the file to be read.")
      .def("setFileName", &VTKReader<T>::setFileName,
           "Set the name of the input file.")
      .def("apply", &VTKReader<T>::apply, "Read the mesh.");

  // VTKWriter
  pybind11::class_<VTKWriter<T>, SmartPointer<VTKWriter<T>>>(module,
                                                             "VTKWriter")
      // constructors
      .def(pybind11::init(&SmartPointer<VTKWriter<T>>::New<>))
      .def(pybind11::init(
          &SmartPointer<VTKWriter<T>>::New<SmartPointer<Mesh<T>> &>))
      .def(pybind11::init(
          &SmartPointer<VTKWriter<T>>::New<SmartPointer<Mesh<T>> &,
                                           std::string>))
      .def(pybind11::init([](SmartPointer<Mesh<T>> &mesh, FileFormatEnum format,
                             std::string s) {
        return SmartPointer<VTKWriter<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &VTKWriter<T>::setMesh, "Set the mesh to output.")
      .def("setFileFormat", &VTKWriter<T>::setFileFormat,
           "Set the file format, the mesh should be written to.")
      .def("setFileName", &VTKWriter<T>::setFileName,
           "Set the name of the output file.")
      .def("apply", &VTKWriter<T>::apply, "Write the mesh.");

  // Writer
  pybind11::class_<Writer<T, D>, SmartPointer<Writer<T, D>>>(module, "Writer")
      // constructors
      .def(pybind11::init(&SmartPointer<Writer<T, D>>::New<>))
      .def(pybind11::init(
          &SmartPointer<Writer<T, D>>::New<SmartPointer<Domain<T, D>> &>))
      .def(pybind11::init(
          &SmartPointer<Writer<T, D>>::New<SmartPointer<Domain<T, D>> &,
                                           std::string>))
      // methods
      .def("setLevelSet", &Writer<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &Writer<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &Writer<T, D>::apply, "Write to file.");

// WriteVisualizationMesh
#ifdef VIENNALS_USE_VTK
  pybind11::class_<WriteVisualizationMesh<T, D>,
                   SmartPointer<WriteVisualizationMesh<T, D>>>(
      module, "WriteVisualizationMesh")
      // constructors
      .def(pybind11::init(&SmartPointer<WriteVisualizationMesh<T, D>>::New<>))
      .def(pybind11::init(&SmartPointer<WriteVisualizationMesh<T, D>>::New<
                          SmartPointer<Domain<T, D>> &>))
      // methods
      .def("insertNextLevelSet",
           &WriteVisualizationMesh<T, D>::insertNextLevelSet,
           "Insert next level set to convert. Bigger level sets wrapping "
           "smaller ones, should be inserted last.")
      .def("setFileName", &WriteVisualizationMesh<T, D>::setFileName,
           "Set Name of File to write.")
      .def("setExtractHullMesh",
           &WriteVisualizationMesh<T, D>::setExtractHullMesh,
           "Whether to extract a hull mesh. Defaults to false.")
      .def("setExtractVolumeMesh",
           &WriteVisualizationMesh<T, D>::setExtractVolumeMesh,
           " Whether to extract a tetra volume mesh. Defaults to true.")
      .def("apply", &WriteVisualizationMesh<T, D>::apply,
           "Make and write mesh.");
#endif

  // 2D Domain in 3D import and 3D domain in 2D import
  //   constexpr int dim = VIENNALS_PYTHON_DIMENSION == 2 ? 3 : 2;
  //   pybind11::class_<Domain<T, dim>, SmartPointer<Domain<T, dim>>>(
  //       module, ("Domain" + std::to_string(dim) + "D").c_str())
  //       // constructors
  //       .def(pybind11::init(&SmartPointer<Domain<T, dim>>::New<>))
  //       .def(pybind11::init(
  //                &SmartPointer<Domain<T, dim>>::New<viennahrle::CoordType>),
  //            pybind11::arg("gridDelta") = 1.0)
  //       .def(pybind11::init([](std::array<viennahrle::CoordType, 2 * D>
  //       bounds,
  //                              std::array<BoundaryConditionEnum, D> bcs,
  //                              viennahrle::CoordType gridDelta) {
  //              return SmartPointer<Domain<T, dim>>::New(bounds.data(),
  //              bcs.data(),
  //                                                       gridDelta);
  //            }),
  //            pybind11::arg("bounds"), pybind11::arg("boundaryConditions"),
  //            pybind11::arg("gridDelta") = 1.0)
  //       .def(pybind11::init(&SmartPointer<Domain<T, dim>>::New<
  //                           std::vector<viennahrle::CoordType>,
  //                           std::vector<unsigned>, viennahrle::CoordType>),
  //            pybind11::arg("bounds"), pybind11::arg("boundaryConditions"),
  //            pybind11::arg("gridDelta") = 1.0)
  //       .def(pybind11::init(
  //           &SmartPointer<Domain<T, dim>>::New<SmartPointer<Domain<T, dim>>
  //           &>))
  //       .def(pybind11::init(
  //           &SmartPointer<Domain<T, dim>>::New<viennahrle::Grid<dim> &>))
  //       // methods
  //       .def("deepCopy", &Domain<T, dim>::deepCopy,
  //            "Copy lsDomain in this lsDomain.")
  //       .def("getNumberOfSegments", &Domain<T, dim>::getNumberOfSegments,
  //            "Get the number of segments, the level set structure is divided
  //            " "into.")
  //       .def("getNumberOfPoints", &Domain<T, dim>::getNumberOfPoints,
  //            "Get the number of defined level set values.")
  //       .def("getLevelSetWidth", &Domain<T, dim>::getLevelSetWidth,
  //            "Get the number of layers of level set points around the
  //            explicit " "surface.")
  //       .def("setLevelSetWidth", &Domain<T, dim>::setLevelSetWidth,
  //            "Set the number of layers of level set points which should be "
  //            "stored around the explicit surface.")
  //       .def("clearMetaData", &Domain<T, dim>::clearMetaData,
  //            "Clear all metadata stored in the level set.")
  //       // allow filehandle to be passed and default to python standard
  //       output .def(
  //           "print",
  //           [](Domain<T, D> &d, pybind11::object fileHandle) {
  //             if (!(pybind11::hasattr(fileHandle, "write") &&
  //                   pybind11::hasattr(fileHandle, "flush"))) {
  //               throw pybind11::type_error(
  //                   "MyClass::read_from_file_like_object(file): incompatible
  //                   " "function argument:  `file` must be a file-like object,
  //                   but "
  //                   "`" +
  //                   (std::string)(pybind11::repr(fileHandle)) + "`
  //                   provided");
  //             }
  //             pybind11::detail::pythonbuf buf(fileHandle);
  //             std::ostream stream(&buf);
  //             d.print(stream);
  //           },
  //           pybind11::arg("stream") =
  //               pybind11::module::import("sys").attr("stdout"));

  // Also wrap hrleGrid so it can be used to create new LevelSets
  pybind11::class_<viennahrle::Grid<D>>(module, "hrleGrid");
}
