#include <pybind11/functional.h>
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
#include <lsCompareChamfer.hpp>
#include <lsCompareCriticalDimensions.hpp>
#include <lsCompareNarrowBand.hpp>
#include <lsCompareSparseField.hpp>
#include <lsCompareVolume.hpp>
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
#include <lsToMultiSurfaceMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsTransformMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKRenderWindow.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>
#include <lsVersion.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <lsWriter.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

using namespace viennals;
namespace py = pybind11;

// always use double for python export
typedef double T;

PYBIND11_DECLARE_HOLDER_TYPE(TemplateType, SmartPointer<TemplateType>);

// GeometricAdvectDistribution
template <int D>
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

template <int D> void bindApi(py::module &module) {
  // Also wrap hrleGrid so it can be used to create new LevelSets
  py::class_<viennahrle::Grid<D>>(module, "hrleGrid");

  // Domain
  py::class_<Domain<T, D>, SmartPointer<Domain<T, D>>>(module, "Domain")
      // constructors
      .def(py::init(&SmartPointer<Domain<T, D>>::template New<>))
      .def(
          py::init(
              &SmartPointer<Domain<T, D>>::template New<viennahrle::CoordType>),
          py::arg("gridDelta") = 1.0)
      //  .def(py::init(
      //      &SmartPointer<Domain<T, D>>::New<viennahrle::CoordType *,
      //                                       BoundaryConditionEnum *>))
      .def(py::init([](std::array<viennahrle::CoordType, 2 * D> bounds,
                       std::array<BoundaryConditionEnum, D> bcs,
                       viennahrle::CoordType gridDelta) {
             return SmartPointer<Domain<T, D>>::New(bounds.data(), bcs.data(),
                                                    gridDelta);
           }),
           py::arg("bounds"), py::arg("boundaryConditions"),
           py::arg("gridDelta") = 1.0)
      .def(py::init(&SmartPointer<Domain<T, D>>::template New<
                    std::vector<viennahrle::CoordType>, std::vector<unsigned>,
                    viennahrle::CoordType>),
           py::arg("bounds"), py::arg("boundaryConditions"),
           py::arg("gridDelta") = 1.0)
      //  .def(py::init(
      //      &SmartPointer<Domain<T, D>>::New<Domain<T,
      //      D>::PointValueVectorType,
      //                                       viennahrle::CoordType *,
      //                                       BoundaryConditionEnum *>))
      //  .def(py::init(&SmartPointer<Domain<T, D>>::New<
      //                      Domain<T, D>::PointValueVectorType,
      //                      viennahrle::CoordType
      //                      *, BoundaryConditionEnum *,
      //                      viennahrle::CoordType>))
      .def(py::init(&SmartPointer<Domain<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(
          &SmartPointer<Domain<T, D>>::template New<viennahrle::Grid<D> &>))
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
          [](Domain<T, D> &d, py::object fileHandle) {
            if (!(py::hasattr(fileHandle, "write") &&
                  py::hasattr(fileHandle, "flush"))) {
              throw py::type_error(
                  "MyClass::read_from_file_like_object(file): incompatible "
                  "function argument:  `file` must be a file-like object, but "
                  "`" +
                  (std::string)(py::repr(fileHandle)) + "` provided");
            }
            py::detail::pythonbuf buf(fileHandle);
            std::ostream stream(&buf);
            d.print(stream);
          },
          py::arg("stream") = py::module::import("sys").attr("stdout"));

  // Advect
  py::class_<Advect<T, D>, SmartPointer<Advect<T, D>>>(module, "Advect")
      // constructors
      .def(py::init(&SmartPointer<Advect<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Advect<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(
          &SmartPointer<Advect<T, D>>::template New<
              SmartPointer<Domain<T, D>> &, SmartPointer<VelocityField<T>> &>))
      // getters and setters
      .def("insertNextLevelSet", &Advect<T, D>::insertNextLevelSet,
           "Insert next level set to use for advection.")
      .def("clearLevelSets", &Advect<T, D>::clearLevelSets,
           "Clear all level sets used for advection.")
      .def("setVelocityField", &Advect<T, D>::setVelocityField,
           "Set the velocity to use for advection.")
      .def("setAdvectionTime", &Advect<T, D>::setAdvectionTime,
           "Set the time until when the level set should be advected.")
      .def("setSingleStep", &Advect<T, D>::setSingleStep, py::arg("singleStep"),
           "Set whether only a single advection step should be performed.")
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
      .def("setAdaptiveTimeStepping", &Advect<T, D>::setAdaptiveTimeStepping,
           py::arg("enabled") = true, py::arg("subdivisions") = 20,
           "Enable/disable adaptive time stepping and set the number of "
           "subdivisions.")
      .def(
          "setSaveAdvectionVelocities",
          &Advect<T, D>::setSaveAdvectionVelocities,
          "Set whether the velocities applied to each point should be saved in "
          "the level set for debug purposes.")
      .def("setCheckDissipation", &Advect<T, D>::setCheckDissipation,
           py::arg("check"), "Enable/disable dissipation checking.")
      .def("setUpdatePointData", &Advect<T, D>::setUpdatePointData,
           py::arg("update"),
           "Enable/disable updating point data after advection.")
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
      .def("setSpatialScheme", &Advect<T, D>::setSpatialScheme,
           "Set the spatial discretization scheme to use during advection.")
      .def("setTemporalScheme", &Advect<T, D>::setTemporalScheme,
           "Set the time integration scheme to use during advection.")
      .def("setIntegrationScheme", &Advect<T, D>::setIntegrationScheme,
           "(DEPRECATED, use setSpatialScheme instead) Set the spatial "
           "discretization scheme to use during advection.")
      .def("setDissipationAlpha", &Advect<T, D>::setDissipationAlpha,
           "Set the dissipation value to use for Lax Friedrichs spatial "
           "discretization.")
      .def("setUpdatePointData", &Advect<T, D>::setUpdatePointData,
           "Set whether the point data in the old LS should be translated to "
           "the advected LS. Defaults to true.")
      .def(
          "setVelocityUpdateCallback", &Advect<T, D>::setVelocityUpdateCallback,
          "Set a callback function that is called after the level set has been "
          "updated during intermediate time integration steps (e.g. RK2, RK3).")
      .def("prepareLS", &Advect<T, D>::prepareLS, "Prepare the level-set.")
      // need scoped release since we are calling a python method from
      // parallelised C++ code here
      .def("apply", &Advect<T, D>::apply,
           py::call_guard<py::gil_scoped_release>(), "Perform advection.");

  py::class_<lsInternal::StencilLocalLaxFriedrichsScalar<T, D, 1>>(
      module, "StencilLocalLaxFriedrichsScalar")
      .def_static(
          "setMaxDissipation",
          &lsInternal::StencilLocalLaxFriedrichsScalar<T, D,
                                                       1>::setMaxDissipation,
          py::arg("maxDissipation"));

  module.def("PrepareStencilLocalLaxFriedrichs",
             &PrepareStencilLocalLaxFriedrichs<T, D>, py::arg("levelSets"),
             py::arg("isDepo"));

  module.def("FinalizeStencilLocalLaxFriedrichs",
             &FinalizeStencilLocalLaxFriedrichs<T, D>, py::arg("levelSets"));

  // BooleanOperation
  py::class_<BooleanOperation<T, D>, SmartPointer<BooleanOperation<T, D>>>(
      module, "BooleanOperation")
      // constructors
      .def(py::init(&SmartPointer<BooleanOperation<T, D>>::template New<>))
      .def(py::init(&SmartPointer<BooleanOperation<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(
          py::init(&SmartPointer<BooleanOperation<T, D>>::template New<
                   SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(py::init(
          [](SmartPointer<Domain<T, D>> &domain, BooleanOperationEnum op) {
            return SmartPointer<BooleanOperation<T, D>>::New(domain, op);
          }))
      .def(py::init([](SmartPointer<Domain<T, D>> &domainA,
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

  py::class_<CalculateCurvatures<T, D>,
             SmartPointer<CalculateCurvatures<T, D>>>(module,
                                                      "CalculateCurvatures")
      // constructors
      .def(py::init(&SmartPointer<CalculateCurvatures<T, D>>::template New<>))
      .def(py::init(&SmartPointer<CalculateCurvatures<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(py::init([](SmartPointer<Domain<T, D>> &domain, CurvatureEnum type) {
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

  // CalculateNormalVectors
  py::class_<CalculateNormalVectors<T, D>,
             SmartPointer<CalculateNormalVectors<T, D>>>(
      module, "CalculateNormalVectors")
      // constructors
      .def(
          py::init(&SmartPointer<CalculateNormalVectors<T, D>>::template New<>))
      .def(py::init(&SmartPointer<CalculateNormalVectors<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &CalculateNormalVectors<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("setMaxValue", &CalculateNormalVectors<T, D>::setMaxValue,
           "Set the maximum value for which normals should be calculated.")
      .def("setMethod", &CalculateNormalVectors<T, D>::setMethod,
           "Set the method to use for normal calculation.")
      .def("apply", &CalculateNormalVectors<T, D>::apply,
           "Perform normal vector calculation.");

  // CalculateVisibilities
  py::class_<CalculateVisibilities<T, D>,
             SmartPointer<CalculateVisibilities<T, D>>>(module,
                                                        "CalculateVisibilities")
      .def(py::init(
          &SmartPointer<CalculateVisibilities<T, D>>::template New<
              SmartPointer<Domain<T, D>> &, const Vec3D<T> &, std::string>))
      .def("apply", &CalculateVisibilities<T, D>::apply);

  // Check
  py::class_<Check<T, D>, SmartPointer<Check<T, D>>>(module, "Check")
      // constructors
      .def(py::init(&SmartPointer<Check<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Check<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &Check<T, D>::setLevelSet,
           "Set levelset for which to calculate normal vectors.")
      .def("apply", &Check<T, D>::apply, "Perform check.");

  // PointCloud
  py::class_<PointCloud<T, D>, SmartPointer<PointCloud<T, D>>>(module,
                                                               "PointCloud")
      // constructors
      .def(py::init(&SmartPointer<PointCloud<T, D>>::template New<
                    const std::vector<VectorType<T, D>> &>))
      // methods
      .def("insertNextPoint",
           (void(PointCloud<T, D>::*)(const VectorType<T, D> &)) &
               PointCloud<T, D>::insertNextPoint);

  // ConvexHull
  py::class_<ConvexHull<T, D>, SmartPointer<ConvexHull<T, D>>>(module,
                                                               "ConvexHull")
      // constructors
      .def(py::init(&SmartPointer<ConvexHull<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ConvexHull<T, D>>::template New<
                    SmartPointer<Mesh<T>> &, SmartPointer<PointCloud<T, D>> &>))
      // methods
      .def("setMesh", &ConvexHull<T, D>::setMesh,
           "Set mesh object where the generated mesh should be stored.")
      .def("setPointCloud", &ConvexHull<T, D>::setPointCloud,
           "Set point cloud used to generate mesh.")
      .def("apply", &ConvexHull<T, D>::apply, "Generate Hull.");

  // DetectFeatures
  py::class_<DetectFeatures<T, D>, SmartPointer<DetectFeatures<T, D>>>(
      module, "DetectFeatures")
      // constructors
      .def(py::init(&SmartPointer<DetectFeatures<T, D>>::template New<>))
      .def(py::init(&SmartPointer<DetectFeatures<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<DetectFeatures<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, T>))
      // some constructors need lambda to work: seems to be an issue with
      // implicit move constructor
      .def(py::init([](SmartPointer<Domain<T, D>> &domain, T maxValue,
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

  // GeometricAdvect
  py::class_<GeometricAdvect<T, D>, SmartPointer<GeometricAdvect<T, D>>>(
      module, "GeometricAdvect")
      // constructors
      .def(py::init(&SmartPointer<GeometricAdvect<T, D>>::template New<>))
      .def(py::init(
          &SmartPointer<GeometricAdvect<T, D>>::template New<
              SmartPointer<Domain<T, D>> &,
              SmartPointer<
                  GeometricAdvectDistribution<viennahrle::CoordType, D>> &>))
      // methods
      .def("setLevelSet", &GeometricAdvect<T, D>::setLevelSet,
           "Set levelset to advect.")
      .def(
          "setAdvectionDistribution",
          &GeometricAdvect<T, D>::setAdvectionDistribution,
          "Set advection distribution to use as kernel for the fast advection.")
      .def("apply", &GeometricAdvect<T, D>::apply,
           py::call_guard<py::gil_scoped_release>(), "Perform advection.");

  // GeometricAdvectDistributions
  py::class_<GeometricAdvectDistribution<T, D>,
             SmartPointer<GeometricAdvectDistribution<T, D>>,
             PylsGeometricAdvectDistribution<D>>(module,
                                                 "GeometricAdvectDistribution")
      // constructors
      .def(py::init<>())
      // methods
      .def("isInside", &GeometricAdvectDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance",
           &GeometricAdvectDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &GeometricAdvectDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.")
      .def("prepare", &GeometricAdvectDistribution<T, D>::prepare,
           "Prepare the distribution for use with the passed level set.")
      .def("finalize", &GeometricAdvectDistribution<T, D>::finalize,
           "Finalize the distribution after use with the level set.");

  py::class_<SphereDistribution<T, D>, SmartPointer<SphereDistribution<T, D>>,
             GeometricAdvectDistribution<T, D>>(module, "SphereDistribution")
      // constructors
      .def(py::init(&SmartPointer<SphereDistribution<T, D>>::template New<T>))
      // methods
      .def("isInside", &SphereDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &SphereDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &SphereDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  py::class_<BoxDistribution<T, D>, SmartPointer<BoxDistribution<T, D>>,
             GeometricAdvectDistribution<T, D>>(module, "BoxDistribution")
      // constructors
      .def(py::init(&SmartPointer<BoxDistribution<T, D>>::template New<
                    const std::array<T, 3>>))
      // methods
      .def("isInside", &BoxDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance", &BoxDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &BoxDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  py::class_<CustomSphereDistribution<T, D>,
             SmartPointer<CustomSphereDistribution<T, D>>,
             GeometricAdvectDistribution<T, D>>(module,
                                                "CustomSphereDistribution")
      // constructors
      .def(py::init(&SmartPointer<CustomSphereDistribution<T, D>>::template New<
                    const std::vector<T> &>))
      // methods
      .def("isInside", &CustomSphereDistribution<T, D>::isInside,
           "Check whether passed point is inside the distribution.")
      .def("getSignedDistance",
           &CustomSphereDistribution<T, D>::getSignedDistance,
           "Get the signed distance of the passed point to the surface of the "
           "distribution.")
      .def("getBounds", &CustomSphereDistribution<T, D>::getBounds,
           "Get the cartesian bounds of the distribution.");

  // Expand
  py::class_<Expand<T, D>, SmartPointer<Expand<T, D>>>(module, "Expand")
      // constructors
      .def(py::init(&SmartPointer<Expand<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Expand<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<Expand<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, int>))
      // methods
      .def("setLevelSet", &Expand<T, D>::setLevelSet, "Set levelset to expand.")
      .def("setWidth", &Expand<T, D>::setWidth, "Set the width to expand to.")
      .def("apply", &Expand<T, D>::apply, "Perform expansion.");

  // FromSurfaceMesh
  py::class_<FromSurfaceMesh<T, D>, SmartPointer<FromSurfaceMesh<T, D>>>(
      module, "FromSurfaceMesh")
      // constructors
      .def(py::init(&SmartPointer<FromSurfaceMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<FromSurfaceMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      .def(py::init(
          &SmartPointer<FromSurfaceMesh<T, D>>::template New<
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
  py::class_<FromVolumeMesh<T, D>, SmartPointer<FromVolumeMesh<T, D>>>(
      module, "FromVolumeMesh")
      // constructors
      .def(py::init(&SmartPointer<FromVolumeMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<FromVolumeMesh<T, D>>::template New<
                    typename FromVolumeMesh<T, D>::GridType &,
                    SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<FromVolumeMesh<T, D>>::template New<
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
  py::class_<FromMesh<T, D>, SmartPointer<FromMesh<T, D>>>(module, "FromMesh")
      .def(py::init(&SmartPointer<FromMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<FromMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      .def("setMesh", &FromMesh<T, D>::setMesh, "Set the mesh to read from.")
      .def("setSortPointList", &FromMesh<T, D>::setSortPointList)
      .def("apply", &FromMesh<T, D>::apply);

  // lsGeometries
  // Sphere
  py::class_<Sphere<T, D>, SmartPointer<Sphere<T, D>>>(module, "Sphere")
      // constructors
      .def(py::init(&SmartPointer<Sphere<T, D>>::template New<
                    const std::vector<T> & /*origin*/, T /*radius*/>),
           py::arg("origin"), py::arg("radius"));
  // Plane
  py::class_<Plane<T, D>, SmartPointer<Plane<T, D>>>(module, "Plane")
      // constructors
      .def(py::init(&SmartPointer<Plane<T, D>>::template New<
                    const std::vector<T> & /*origin*/,
                    const std::vector<T> & /*normal*/>),
           py::arg("origin"), py::arg("normal"));
  // Box
  py::class_<Box<T, D>, SmartPointer<Box<T, D>>>(module, "Box")
      // constructors
      .def(py::init(&SmartPointer<Box<T, D>>::template New<
                    const std::vector<T> & /*minPoint*/,
                    const std::vector<T> & /*maxPoint*/>),
           py::arg("minPoint"), py::arg("maxPoint"));
  // Cylinder
  py::class_<Cylinder<T, D>, SmartPointer<Cylinder<T, D>>>(module, "Cylinder")
      // constructors
      .def(
          py::init(&SmartPointer<Cylinder<T, D>>::template New<
                   const std::vector<T> & /*origin*/,
                   const std::vector<T> & /*axisDirection*/, const T /*height*/,
                   const T /*radius*/, const T /*topRadius*/>),
          py::arg("origin"), py::arg("axisDirection"), py::arg("height"),
          py::arg("radius"), py::arg("topRadius") = 0.);

  // MakeGeometry
  py::class_<MakeGeometry<T, D>, SmartPointer<MakeGeometry<T, D>>>(
      module, "MakeGeometry")
      // constructors
      .def(py::init(&SmartPointer<MakeGeometry<T, D>>::template New<>))
      .def(py::init(&SmartPointer<MakeGeometry<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(
          py::init(&SmartPointer<MakeGeometry<T, D>>::template New<
                   SmartPointer<Domain<T, D>> &, SmartPointer<Sphere<T, D>> &>))
      .def(py::init(&SmartPointer<MakeGeometry<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Plane<T, D>> &>))
      .def(py::init(&SmartPointer<MakeGeometry<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Box<T, D>> &>))
      .def(py::init(
          &SmartPointer<MakeGeometry<T, D>>::template New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Cylinder<T, D>> &>))
      .def(py::init(
          &SmartPointer<MakeGeometry<T, D>>::template New<
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
  py::class_<MarkVoidPoints<T, D>, SmartPointer<MarkVoidPoints<T, D>>>(
      module, "MarkVoidPoints")
      // constructors
      .def(py::init(&SmartPointer<MarkVoidPoints<T, D>>::template New<>))
      .def(py::init(&SmartPointer<MarkVoidPoints<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<MarkVoidPoints<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, bool &>))
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
      .def("setSaveComponentIds", &MarkVoidPoints<T, D>::setSaveComponentIds,
           "Save the connectivity information of all LS points in the "
           "pointData of the level set.")
      .def("getNumberOfComponents",
           &MarkVoidPoints<T, D>::getNumberOfComponents,
           "Get the number of connected components found in the level set.")
      .def("apply", &MarkVoidPoints<T, D>::apply, "Mark void points.");

  // Prune
  py::class_<Prune<T, D>, SmartPointer<Prune<T, D>>>(module, "Prune")
      // constructors
      .def(py::init(&SmartPointer<Prune<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Prune<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &Prune<T, D>::setLevelSet, "Set levelset to prune.")
      .def("apply", &Prune<T, D>::apply, "Perform pruning operation.");

  // Reader
  py::class_<Reader<T, D>, SmartPointer<Reader<T, D>>>(module, "Reader")
      // constructors
      .def(py::init(&SmartPointer<Reader<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Reader<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<Reader<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, std::string>))
      // methods
      .def("setLevelSet", &Reader<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &Reader<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &Reader<T, D>::apply, "Write to file.");

  // Reduce
  py::class_<Reduce<T, D>, SmartPointer<Reduce<T, D>>>(module, "Reduce")
      // constructors
      .def(py::init(&SmartPointer<Reduce<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Reduce<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<Reduce<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, int>))
      .def(py::init(&SmartPointer<Reduce<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, int, bool>))
      // methods
      .def("setLevelSet", &Reduce<T, D>::setLevelSet, "Set levelset to reduce.")
      .def("setWidth", &Reduce<T, D>::setWidth, "Set the width to reduce to.")
      .def("setNoNewSegment", &Reduce<T, D>::setNoNewSegment,
           "Set whether the levelset should be segmented anew (balanced across "
           "cores) after reduction.")
      .def("apply", &Reduce<T, D>::apply, "Perform reduction.");

  // RemoveStrayPoints
  py::class_<RemoveStrayPoints<T, D>, SmartPointer<RemoveStrayPoints<T, D>>>(
      module, "RemoveStrayPoints")
      // constructors
      .def(py::init(&SmartPointer<RemoveStrayPoints<T, D>>::template New<>))
      .def(py::init(&SmartPointer<RemoveStrayPoints<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSet", &RemoveStrayPoints<T, D>::setLevelSet,
           "Set levelset for stray point removal.")
      .def("setVoidTopSurface", &RemoveStrayPoints<T, D>::setVoidTopSurface,
           "Set the logic by which to choose the surface which should be kept. "
           "All other LS values will be marked as stray points and removed.")
      .def("apply", &RemoveStrayPoints<T, D>::apply, "Remove stray points.");

  // ToDiskMesh
  py::class_<ToDiskMesh<T, D>, SmartPointer<ToDiskMesh<T, D>>>(module,
                                                               "ToDiskMesh")
      // constructors
      .def(py::init(&SmartPointer<ToDiskMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ToDiskMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      // methods
      .def("setLevelSet", &ToDiskMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("clearLevelSets", &ToDiskMesh<T, D>::clearLevelSets,
           "Clear all inserted level sets.")
      .def("insertNextLevelSet", &ToDiskMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the disk mesh.")
      .def("setMesh", &ToDiskMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("setMaterialMap", &ToDiskMesh<T, D>::setMaterialMap,
           "Set the material map to use for the disk mesh.")
      .def("setMaxValue", &ToDiskMesh<T, D>::setMaxValue,
           "Set the maximum level set value to include in the disk mesh.")
      .def("apply", &ToDiskMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToMesh
  py::class_<ToMesh<T, D>, SmartPointer<ToMesh<T, D>>>(module, "ToMesh")
      // constructors
      .def(py::init(&SmartPointer<ToMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ToMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      .def(py::init(
          &SmartPointer<ToMesh<T, D>>::template New<
              SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &, bool>))
      .def(py::init(&SmartPointer<ToMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &, bool,
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
  py::class_<ToSurfaceMesh<T, D>, SmartPointer<ToSurfaceMesh<T, D>>>(
      module, "ToSurfaceMesh")
      // constructors
      .def(py::init(&SmartPointer<ToSurfaceMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ToSurfaceMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      // methods
      .def("setLevelSet", &ToSurfaceMesh<T, D>::setLevelSet,
           "Set levelset to mesh.")
      .def("setMesh", &ToSurfaceMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("setUpdatePointData", &ToSurfaceMesh<T, D>::setUpdatePointData,
           "Set whether to update point data. Defaults to true.")
      .def("setSharpCorners", &ToSurfaceMesh<T, D>::setSharpCorners,
           "Set whether to preserve sharp corners. Defaults to false.")
      .def("apply", &ToSurfaceMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToMultiSurfaceMesh
  py::class_<ToMultiSurfaceMesh<T, D>, SmartPointer<ToMultiSurfaceMesh<T, D>>>(
      module, "ToMultiSurfaceMesh")
      // constructors
      .def(py::init(
               &SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<double,
                                                                     double>),
           py::arg("eps") = 1e-12, py::arg("minNodeDistFactor") = 0.02)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &,
                    double, double>),
           py::arg("domain"), py::arg("mesh"), py::arg("eps") = 1e-12,
           py::arg("minNodeDistFactor") = 0.02)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    std::vector<SmartPointer<Domain<T, D>>> &,
                    SmartPointer<Mesh<T>> &, double, double>),
           py::arg("domains"), py::arg("mesh"), py::arg("eps") = 1e-12,
           py::arg("minNodeDistFactor") = 0.02)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &, double, double>),
           py::arg("mesh"), py::arg("eps") = 1e-12,
           py::arg("minNodeDistFactor") = 0.02)
      // methods
      .def("insertNextLevelSet", &ToMultiSurfaceMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the mesh.")
      .def("clearLevelSets", &ToMultiSurfaceMesh<T, D>::clearLevelSets,
           "Clear all inserted level sets.")
      .def("setMesh", &ToMultiSurfaceMesh<T, D>::setMesh,
           "Set the mesh to generate.")
      .def("setMaterialMap", &ToMultiSurfaceMesh<T, D>::setMaterialMap,
           "Set the material map to use for the multi surface mesh.")
      .def("setSharpCorners", &ToMultiSurfaceMesh<T, D>::setSharpCorners,
           "Set whether to preserve sharp corners. Defaults to false.")
      .def("apply", &ToMultiSurfaceMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // ToVoxelMesh
  py::class_<ToVoxelMesh<T, D>, SmartPointer<ToVoxelMesh<T, D>>>(module,
                                                                 "ToVoxelMesh")
      // constructors
      .def(py::init(&SmartPointer<ToVoxelMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ToVoxelMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<ToVoxelMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<ToVoxelMesh<T, D>>::template New<
                    std::vector<SmartPointer<Domain<T, D>>> &,
                    SmartPointer<Mesh<T>> &>))
      // methods
      .def("insertNextLevelSet", &ToVoxelMesh<T, D>::insertNextLevelSet,
           "Insert next level set to output in the mesh.")
      .def("clearLevelSets", &ToVoxelMesh<T, D>::clearLevelSets,
           "Clear all inserted level sets.")
      .def("setMesh", &ToVoxelMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("apply", &ToVoxelMesh<T, D>::apply,
           "Convert the levelset to a surface mesh.");

  // Writer
  py::class_<Writer<T, D>, SmartPointer<Writer<T, D>>>(module, "Writer")
      // constructors
      .def(py::init(&SmartPointer<Writer<T, D>>::template New<>))
      .def(py::init(&SmartPointer<Writer<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &>))
      .def(py::init(&SmartPointer<Writer<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, std::string>))
      // methods
      .def("setLevelSet", &Writer<T, D>::setLevelSet,
           "Set levelset to write to file.")
      .def("setFileName", &Writer<T, D>::setFileName,
           "Set the filename for the output file.")
      .def("apply", &Writer<T, D>::apply, "Write to file.");

  // CompareSparseField
  py::class_<CompareSparseField<T, D>, SmartPointer<CompareSparseField<T, D>>>(
      module, "CompareSparseField")
      // constructors
      .def(py::init(&SmartPointer<CompareSparseField<T, D>>::template New<>))
      .def(
          py::init(&SmartPointer<CompareSparseField<T, D>>::template New<
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
      .def("setZRange", &CompareSparseField<T, D>::setZRange,
           "Set the z-coordinate range to restrict the comparison area")
      .def("clearZRange", &CompareSparseField<T, D>::clearZRange,
           "Clear the z-range restriction")
      .def("setOutputMesh", &CompareSparseField<T, D>::setOutputMesh,
           "Set the output mesh where difference values will be stored")
      .def("setFillIteratedWithDistances",
           &CompareSparseField<T, D>::setFillIteratedWithDistances,
           "Set whether to fill the iterated level set with distance values")
      .def("setExpandedLevelSetWidth",
           &CompareSparseField<T, D>::setExpandedLevelSetWidth,
           "Set the expansion width for the expanded level set")
      .def("apply", &CompareSparseField<T, D>::apply,
           "Apply the comparison and calculate the sum of squared "
           "differences.")
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

// WriteVisualizationMesh
#ifdef VIENNALS_USE_VTK
  py::class_<WriteVisualizationMesh<T, D>,
             SmartPointer<WriteVisualizationMesh<T, D>>>(
      module, "WriteVisualizationMesh")
      // constructors
      .def(
          py::init(&SmartPointer<WriteVisualizationMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<WriteVisualizationMesh<T, D>>::template New<
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
      .def("setMetaData", &WriteVisualizationMesh<T, D>::setMetaData,
           "Set the metadata to be written to the file.")
      .def("addMetaData",
           py::overload_cast<const std::string &, T>(
               &WriteVisualizationMesh<T, D>::addMetaData),
           "Add a single metadata entry to the file.")
      .def("addMetaData",
           py::overload_cast<const std::string &, const std::vector<T> &>(
               &WriteVisualizationMesh<T, D>::addMetaData),
           "Add a single metadata entry to the file.")
      .def("addMetaData",
           py::overload_cast<
               const std::unordered_map<std::string, std::vector<T>> &>(
               &WriteVisualizationMesh<T, D>::addMetaData),
           "Add metadata to the file.")
      .def("apply", &WriteVisualizationMesh<T, D>::apply,
           "Make and write mesh.");
#endif

  // CompareVolume
  py::class_<CompareVolume<T, D>, SmartPointer<CompareVolume<T, D>>>(
      module, "CompareVolume")
      // constructors
      .def(py::init(&SmartPointer<CompareVolume<T, D>>::template New<>))
      .def(
          py::init(&SmartPointer<CompareVolume<T, D>>::template New<
                   SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetTarget", &CompareVolume<T, D>::setLevelSetTarget,
           "Sets the target level set.")
      .def("setLevelSetSample", &CompareVolume<T, D>::setLevelSetSample,
           "Sets the sample level set.")
      .def("setDefaultIncrement", &CompareVolume<T, D>::setDefaultIncrement,
           "Set default increment value")
      .def("setXRangeAndIncrement", &CompareVolume<T, D>::setXRangeAndIncrement,
           "Sets the x-range and custom increment value")
      .def("setYRangeAndIncrement", &CompareVolume<T, D>::setYRangeAndIncrement,
           "Sets the y-range and custom increment value")
      .def("setZRangeAndIncrement", &CompareVolume<T, D>::setZRangeAndIncrement,
           "Sets the z-range and custom increment value")
      .def("setOutputMesh", &CompareVolume<T, D>::setOutputMesh,
           "Set the output mesh where difference areas will be stored")
      .def("getVolumeMismatch", &CompareVolume<T, D>::getVolumeMismatch,
           "Returns the computed volume mismatch.")
      .def("getAreaMismatch", &CompareVolume<T, D>::getAreaMismatch,
           "Returns the computed area mismatch.")
      .def("getCustomVolumeMismatch",
           &CompareVolume<T, D>::getCustomVolumeMismatch,
           "Returns the computed volume mismatch, with custom increments "
           "applied.")
      .def("getCustomAreaMismatch", &CompareVolume<T, D>::getCustomAreaMismatch,
           "Returns the computed area mismatch, with custom increments "
           "applied.")
      .def("getCellCount", &CompareVolume<T, D>::getCellCount,
           "Returns the number of cells where the level sets differ.")
      .def("getCustomCellCount", &CompareVolume<T, D>::getCustomCellCount,
           "Returns the number of cells where the level sets differ, with "
           "custom increments applied.")
      .def("apply", &CompareVolume<T, D>::apply,
           "Computes the volume difference between the two level sets.");

  if constexpr (D == 2) {
    module.attr("CompareArea") = module.attr("CompareVolume");
  }

  // CompareChamfer
  py::class_<CompareChamfer<T, D>, SmartPointer<CompareChamfer<T, D>>>(
      module, "CompareChamfer")
      // constructors
      .def(py::init(&SmartPointer<CompareChamfer<T, D>>::template New<>))
      .def(
          py::init(&SmartPointer<CompareChamfer<T, D>>::template New<
                   SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
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
           "Get the forward distance (average distance from target to "
           "sample).")
      .def("getBackwardDistance", &CompareChamfer<T, D>::getBackwardDistance,
           "Get the backward distance (average distance from sample to "
           "target).")
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

  // CompareCriticalDimensions
  py::class_<CompareCriticalDimensions<T, D>,
             SmartPointer<CompareCriticalDimensions<T, D>>>(
      module, "CompareCriticalDimensions")
      // constructors
      .def(py::init(
          &SmartPointer<CompareCriticalDimensions<T, D>>::template New<>))
      .def(
          py::init(&SmartPointer<CompareCriticalDimensions<T, D>>::template New<
                   SmartPointer<Domain<T, D>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setLevelSetTarget",
           &CompareCriticalDimensions<T, D>::setLevelSetTarget,
           "Sets the target level set.")
      .def("setLevelSetSample",
           &CompareCriticalDimensions<T, D>::setLevelSetSample,
           "Sets the sample level set.")
      .def("addRange", &CompareCriticalDimensions<T, D>::addRange,
           py::arg("measureDimension"), py::arg("minBounds"),
           py::arg("maxBounds"), py::arg("findMaximum") = true,
           "Add a generic range specification.")
      .def("addXRange", &CompareCriticalDimensions<T, D>::addXRange,
           py::arg("minX"), py::arg("maxX"), py::arg("findMaximum") = true,
           "Add an X range to find maximum or minimum Y position.")
      .def("addYRange", &CompareCriticalDimensions<T, D>::addYRange,
           py::arg("minY"), py::arg("maxY"), py::arg("findMaximum") = true,
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
      .def(
          "getCriticalDimensionResult",
          [](CompareCriticalDimensions<T, D> &self, size_t index) {
            T posRef, posCmp, diff;
            bool valid =
                self.getCriticalDimensionResult(index, posRef, posCmp, diff);
            if (valid) {
              return py::make_tuple(true, posRef, posCmp, diff);
            } else {
              return py::make_tuple(false, 0.0, 0.0, 0.0);
            }
          },
          py::arg("index"),
          "Get a specific critical dimension result. Returns (valid, "
          "positionTarget, positionSample, difference).")
      .def("getMeanDifference",
           &CompareCriticalDimensions<T, D>::getMeanDifference,
           "Get mean absolute difference across all valid critical "
           "dimensions.")
      .def("getMaxDifference",
           &CompareCriticalDimensions<T, D>::getMaxDifference,
           "Get maximum difference across all valid critical dimensions.")
      .def("getRMSE", &CompareCriticalDimensions<T, D>::getRMSE,
           "Get RMSE across all valid critical dimensions.")
      .def("getAllDifferences",
           &CompareCriticalDimensions<T, D>::getAllDifferences,
           "Get all valid differences as a list.");

  // CompareNarrowBand
  py::class_<CompareNarrowBand<T, D>, SmartPointer<CompareNarrowBand<T, D>>>(
      module, "CompareNarrowBand")
      // constructors
      .def(py::init(&SmartPointer<CompareNarrowBand<T, D>>::template New<>))
      .def(
          py::init(&SmartPointer<CompareNarrowBand<T, D>>::template New<
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
      .def("setZRange", &CompareNarrowBand<T, D>::setZRange,
           "Set the z-coordinate range to restrict the comparison area")
      .def("clearZRange", &CompareNarrowBand<T, D>::clearZRange,
           "Clear the z-range restriction")
      .def("setOutputMesh", &CompareNarrowBand<T, D>::setOutputMesh,
           "Set the output mesh where difference values will be stored")
      .def("setOutputMeshSquaredDifferences",
           &CompareNarrowBand<T, D>::setOutputMeshSquaredDifferences,
           "Set whether to output squared differences (true) or absolute "
           "differences (false)")
      .def("apply", &CompareNarrowBand<T, D>::apply,
           "Apply the comparison and calculate the sum of squared "
           "differences.")
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
}
