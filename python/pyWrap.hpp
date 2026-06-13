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
#include <lsOxidation.hpp>
#include <lsOxidationModel.hpp>
#include <lsOxidationPresets.hpp>
#include <lsPointData.hpp>
#include <lsPrune.hpp>
#include <lsReader.hpp>
#include <lsReduce.hpp>
#include <lsRemoveStrayPoints.hpp>
#include <lsSlice.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToHullMesh.hpp>
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

inline void bindOxidationSharedTypes(py::module &module) {
  py::class_<OxidationParameters<T>>(module, "OxidationParameters",
                                     py::module_local())
      .def(py::init<>())
      .def_readwrite("diffusionCoefficient",
                     &OxidationParameters<T>::diffusionCoefficient)
      .def_readwrite("reactionRate", &OxidationParameters<T>::reactionRate)
      .def_readwrite("transferCoefficient",
                     &OxidationParameters<T>::transferCoefficient)
      .def_readwrite("equilibriumConcentration",
                     &OxidationParameters<T>::equilibriumConcentration)
      .def_readwrite("oxidantMoleculeDensity",
                     &OxidationParameters<T>::oxidantMoleculeDensity)
      .def_readwrite("expansionCoefficient",
                     &OxidationParameters<T>::expansionCoefficient)
      .def_readwrite("velocitySign", &OxidationParameters<T>::velocitySign)
      .def_readwrite("temperature", &OxidationParameters<T>::temperature)
      .def_readwrite("reactionActivationVolume",
                     &OxidationParameters<T>::reactionActivationVolume)
      .def_readwrite("referencePressure",
                     &OxidationParameters<T>::referencePressure)
      .def_readwrite("diffusionActivationVolume",
                     &OxidationParameters<T>::diffusionActivationVolume)
      .def_readwrite("reactionRateRatio111",
                     &OxidationParameters<T>::reactionRateRatio111)
      .def_readwrite("crystalAxis", &OxidationParameters<T>::crystalAxis)
      .def_readwrite("maskTransferCoefficient",
                     &OxidationParameters<T>::maskTransferCoefficient)
      .def_readwrite("maskConcentration",
                     &OxidationParameters<T>::maskConcentration)
      .def_readwrite("minBoundaryDistance",
                     &OxidationParameters<T>::minBoundaryDistance)
      .def_readwrite("maxIterations", &OxidationParameters<T>::maxIterations)
      .def_readwrite("tolerance", &OxidationParameters<T>::tolerance)
      .def_readwrite("relaxation", &OxidationParameters<T>::relaxation)
      .def_readwrite("maxGridPoints", &OxidationParameters<T>::maxGridPoints)
      .def_readwrite("material", &OxidationParameters<T>::material);

  py::class_<OxidationPresets<T>>(module, "OxidationPresets",
                                  py::module_local())
      .def(py::init<>())
      .def_static("wet1000CDealGrove100",
                  &OxidationPresets<T>::wet1000CDealGrove100)
      .def_static("oxideMechanics1000C",
                  &OxidationPresets<T>::oxideMechanics1000C)
      .def_static("siliconNitrideMask1000C",
                  &OxidationPresets<T>::siliconNitrideMask1000C);

  py::class_<OxidationDeformationParameters<T>>(
      module, "OxidationDeformationParameters", py::module_local())
      .def(py::init<>())
      .def_readwrite("viscosity", &OxidationDeformationParameters<T>::viscosity)
      .def_readwrite("bulkModulus",
                     &OxidationDeformationParameters<T>::bulkModulus)
      .def_readwrite("ambientPressure",
                     &OxidationDeformationParameters<T>::ambientPressure)
      .def_readwrite("pressureRelaxation",
                     &OxidationDeformationParameters<T>::pressureRelaxation)
      .def_readwrite("pressureTolerance",
                     &OxidationDeformationParameters<T>::pressureTolerance)
      .def_readwrite(
          "minMechanicsBoundaryDistance",
          &OxidationDeformationParameters<T>::minMechanicsBoundaryDistance)
      .def_readwrite("shearModulus",
                     &OxidationDeformationParameters<T>::shearModulus)
      .def_readwrite("stressRelaxationTime",
                     &OxidationDeformationParameters<T>::stressRelaxationTime)
      .def_readwrite("stressTimeStep",
                     &OxidationDeformationParameters<T>::stressTimeStep)
      .def_readwrite("harmonicIterations",
                     &OxidationDeformationParameters<T>::harmonicIterations)
      .def_readwrite("mechanicsIterations",
                     &OxidationDeformationParameters<T>::mechanicsIterations)
      .def_readwrite("pressureIterations",
                     &OxidationDeformationParameters<T>::pressureIterations)
      .def_readwrite("stokesIterations",
                     &OxidationDeformationParameters<T>::stokesIterations)
      .def_readwrite("mechanicsTolerance",
                     &OxidationDeformationParameters<T>::mechanicsTolerance)
      .def_readwrite("stokesTolerance",
                     &OxidationDeformationParameters<T>::stokesTolerance)
      .def_readwrite("tolerance", &OxidationDeformationParameters<T>::tolerance)
      .def_readwrite("relaxation",
                     &OxidationDeformationParameters<T>::relaxation)
      .def_readwrite("maxGridPoints",
                     &OxidationDeformationParameters<T>::maxGridPoints)
      .def_readwrite("material", &OxidationDeformationParameters<T>::material);

  py::class_<OxidationMaskParameters<T>>(module, "OxidationMaskParameters",
                                         py::module_local())
      .def(py::init<>())
      .def_readwrite("contactMode", &OxidationMaskParameters<T>::contactMode)
      .def_readwrite("temperature", &OxidationMaskParameters<T>::temperature)
      .def_readwrite("referenceTemperature",
                     &OxidationMaskParameters<T>::referenceTemperature)
      .def_readwrite("referenceViscosity",
                     &OxidationMaskParameters<T>::referenceViscosity)
      .def_readwrite("creepActivationEnergy",
                     &OxidationMaskParameters<T>::creepActivationEnergy)
      .def_readwrite("poissonRatio", &OxidationMaskParameters<T>::poissonRatio)
      .def_readwrite("youngModulus", &OxidationMaskParameters<T>::youngModulus)
      .def_readwrite("unilateralContact",
                     &OxidationMaskParameters<T>::unilateralContact)
      .def_readwrite("relaxation", &OxidationMaskParameters<T>::relaxation)
      .def_readwrite("contactLoadRelaxation",
                     &OxidationMaskParameters<T>::contactLoadRelaxation)
      .def_readwrite("contactReleaseFraction",
                     &OxidationMaskParameters<T>::contactReleaseFraction)
      .def_readwrite("multigridSmootherOmega",
                     &OxidationMaskParameters<T>::multigridSmootherOmega)
      .def_readwrite("tolerance", &OxidationMaskParameters<T>::tolerance)
      .def_readwrite("minBoundaryDistance",
                     &OxidationMaskParameters<T>::minBoundaryDistance)
      .def_readwrite("maxIterations",
                     &OxidationMaskParameters<T>::maxIterations)
      .def_readwrite("maxGridPoints",
                     &OxidationMaskParameters<T>::maxGridPoints)
      .def_readwrite("material", &OxidationMaskParameters<T>::material)
      .def_readwrite("anchorBoundaryDirection",
                     &OxidationMaskParameters<T>::anchorBoundaryDirection)
      .def_readwrite("anchorBoundarySide",
                     &OxidationMaskParameters<T>::anchorBoundarySide)
      .def_readwrite("anchorBoundaryLayers",
                     &OxidationMaskParameters<T>::anchorBoundaryLayers);

  py::class_<OxidationCouplingParameters<T>>(
      module, "OxidationCouplingParameters", py::module_local())
      .def(py::init<>())
      .def_readwrite("maxIterations",
                     &OxidationCouplingParameters<T>::maxIterations)
      .def_readwrite("tolerance", &OxidationCouplingParameters<T>::tolerance)
      .def_readwrite("relaxation", &OxidationCouplingParameters<T>::relaxation);
}

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

  // ReactionBoundarySample — sub-grid crossing info returned by
  // OxidationDiffusion::getReactionBoundarySample
  {
    using RBS = typename OxidationDiffusion<T, D>::ReactionBoundarySample;
    py::class_<RBS>(module, "ReactionBoundarySample")
        .def(py::init<>())
        .def_readwrite("found", &RBS::found)
        .def_readwrite("distance", &RBS::distance)
        .def_readwrite("concentration", &RBS::concentration)
        .def_readwrite("crossingAxis", &RBS::crossingAxis)
        .def_readwrite("crossingOffset", &RBS::crossingOffset)
        .def_property(
            "nodeIndex",
            [](const RBS &s) {
              std::array<viennahrle::IndexType, D> arr{};
              for (unsigned i = 0; i < D; ++i)
                arr[i] = s.nodeIndex[i];
              return arr;
            },
            [](RBS &s, const std::array<viennahrle::IndexType, D> &arr) {
              for (unsigned i = 0; i < D; ++i)
                s.nodeIndex[i] = arr[i];
            });
  }

  py::class_<OxidationDiffusion<T, D>, VelocityField<T>,
             SmartPointer<OxidationDiffusion<T, D>>>(module,
                                                     "OxidationDiffusion")
      .def(py::init([](SmartPointer<Domain<T, D>> &reactionInterface,
                       SmartPointer<Domain<T, D>> &ambientInterface,
                       OxidationParameters<T> parameters) {
             return OxidationDiffusion<T, D>::New(reactionInterface,
                                                  ambientInterface, parameters);
           }),
           py::arg("reactionInterface"), py::arg("ambientInterface"),
           py::arg("parameters") = OxidationParameters<T>())
      .def("setReactionInterface",
           &OxidationDiffusion<T, D>::setReactionInterface)
      .def("setAmbientInterface",
           &OxidationDiffusion<T, D>::setAmbientInterface)
      .def("setMaskInterface", &OxidationDiffusion<T, D>::setMaskInterface,
           py::arg("maskInterface"), py::arg("maskSign") = 1)
      .def("clearMaskInterface", &OxidationDiffusion<T, D>::clearMaskInterface)
      .def("setParameters", &OxidationDiffusion<T, D>::setParameters)
      .def("getParameters", &OxidationDiffusion<T, D>::getParameters)
      .def("setOxideSigns", &OxidationDiffusion<T, D>::setOxideSigns)
      .def(
          "setSolveBounds",
          [](OxidationDiffusion<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setSolveBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def("clearSolveBounds", &OxidationDiffusion<T, D>::clearSolveBounds)
      .def("apply", &OxidationDiffusion<T, D>::apply)
      .def("getConcentration",
           (T(OxidationDiffusion<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDiffusion<T, D>::getConcentration)
      .def("getEffectiveReactionRate",
           (T(OxidationDiffusion<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDiffusion<T, D>::getEffectiveReactionRate)
      .def("clearPressureField", &OxidationDiffusion<T, D>::clearPressureField)
      .def("setPressure",
           (void(OxidationDiffusion<T, D>::*)(const Vec3D<T> &, T)) &
               OxidationDiffusion<T, D>::setPressure)
      .def("getIterations", &OxidationDiffusion<T, D>::getIterations)
      .def("getResidual", &OxidationDiffusion<T, D>::getResidual)
      .def("getNumberOfSolutionNodes",
           &OxidationDiffusion<T, D>::getNumberOfSolutionNodes)
      .def("getReactionBoundarySample",
           &OxidationDiffusion<T, D>::getReactionBoundarySample)
      .def("getScalarVelocityFromSample",
           &OxidationDiffusion<T, D>::getScalarVelocityFromSample)
      .def("markSolved", &OxidationDiffusion<T, D>::markSolved)
      .def("markGeometryChanged",
           &OxidationDiffusion<T, D>::markGeometryChanged)
      .def("writePersistentFields",
           &OxidationDiffusion<T, D>::writePersistentFields)
      .def("getConcentrationCache",
           &OxidationDiffusion<T, D>::getConcentrationCache,
           py::return_value_policy::copy)
      .def("setConcentrationCache",
           &OxidationDiffusion<T, D>::setConcentrationCache);

  py::class_<OxidationDeformation<T, D>, VelocityField<T>,
             SmartPointer<OxidationDeformation<T, D>>>(module,
                                                       "OxidationDeformation")
      .def(
          py::init([](SmartPointer<Domain<T, D>> &reactionInterface,
                      SmartPointer<Domain<T, D>> &ambientInterface,
                      SmartPointer<OxidationDiffusion<T, D>> &diffusionField,
                      OxidationParameters<T> oxidationParameters,
                      OxidationDeformationParameters<T> deformationParameters) {
            return OxidationDeformation<T, D>::New(
                reactionInterface, ambientInterface, diffusionField,
                oxidationParameters, deformationParameters);
          }),
          py::arg("reactionInterface"), py::arg("ambientInterface"),
          py::arg("diffusionField"), py::arg("oxidationParameters"),
          py::arg("deformationParameters") =
              OxidationDeformationParameters<T>())
      .def("setMaskInterface", &OxidationDeformation<T, D>::setMaskInterface,
           py::arg("maskInterface"), py::arg("maskSign") = 1)
      .def("clearMaskInterface",
           &OxidationDeformation<T, D>::clearMaskInterface)
      .def(
          "setSolveBounds",
          [](OxidationDeformation<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setSolveBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def("clearSolveBounds", &OxidationDeformation<T, D>::clearSolveBounds)
      .def("apply", &OxidationDeformation<T, D>::apply)
      .def("getVelocity",
           (Vec3D<T>(OxidationDeformation<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDeformation<T, D>::getVelocity)
      .def("getPressure",
           (T(OxidationDeformation<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDeformation<T, D>::getPressure)
      .def("getStrainTrace",
           (T(OxidationDeformation<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDeformation<T, D>::getStrainTrace)
      .def("getStrainRateTensor",
           (std::array<T, 9>(OxidationDeformation<T, D>::*)(const Vec3D<T> &)
                const) &
               OxidationDeformation<T, D>::getStrainRateTensor)
      .def("getStressTensor", (std::array<T, 9>(OxidationDeformation<T, D>::*)(
                                  const Vec3D<T> &) const) &
                                  OxidationDeformation<T, D>::getStressTensor)
      .def("getVonMisesStress",
           (T(OxidationDeformation<T, D>::*)(const Vec3D<T> &) const) &
               OxidationDeformation<T, D>::getVonMisesStress)
      .def("getIterations", &OxidationDeformation<T, D>::getIterations)
      .def("getResidual", &OxidationDeformation<T, D>::getResidual)
      .def("getNumberOfSolutionNodes",
           &OxidationDeformation<T, D>::getNumberOfSolutionNodes)
      .def("avgExpansionSpeed", &OxidationDeformation<T, D>::avgExpansionSpeed)
      .def("setAmbientInterface",
           &OxidationDeformation<T, D>::setAmbientInterface)
      .def("setOxidationParameters",
           &OxidationDeformation<T, D>::setOxidationParameters)
      .def("setDeformationParameters",
           &OxidationDeformation<T, D>::setDeformationParameters)
      .def("setOxideSigns", &OxidationDeformation<T, D>::setOxideSigns)
      .def("setReactionInterface",
           &OxidationDeformation<T, D>::setReactionInterface)
      .def("setDiffusionField", &OxidationDeformation<T, D>::setDiffusionField)
      .def("setMaskVelocityField",
           &OxidationDeformation<T, D>::setMaskVelocityField)
      .def("clearMaskVelocityField",
           &OxidationDeformation<T, D>::clearMaskVelocityField)
      .def("markGeometryChanged",
           &OxidationDeformation<T, D>::markGeometryChanged)
      .def("writeFieldsToLevelSet",
           &OxidationDeformation<T, D>::writeFieldsToLevelSet);

  py::class_<OxidationMaskBending<T, D>, VelocityField<T>,
             SmartPointer<OxidationMaskBending<T, D>>>(module,
                                                       "OxidationMaskBending")
      .def(py::init(
               [](SmartPointer<OxidationDeformation<T, D>> &deformationField,
                  OxidationMaskParameters<T> maskParameters) {
                 return OxidationMaskBending<T, D>::New(deformationField,
                                                        maskParameters);
               }),
           py::arg("deformationField"),
           py::arg("maskParameters") = OxidationMaskParameters<T>())
      .def(py::init(
               [](SmartPointer<OxidationDeformation<T, D>> &deformationField,
                  SmartPointer<Domain<T, D>> &maskInterface,
                  OxidationMaskParameters<T> maskParameters, int maskSign) {
                 return OxidationMaskBending<T, D>::New(
                     deformationField, maskInterface, maskParameters, maskSign);
               }),
           py::arg("deformationField"), py::arg("maskInterface"),
           py::arg("maskParameters") = OxidationMaskParameters<T>(),
           py::arg("maskSign") = 1)
      .def("setMaskInterface", &OxidationMaskBending<T, D>::setMaskInterface,
           py::arg("maskInterface"), py::arg("maskSign") = 1)
      .def("setAmbientInterface",
           &OxidationMaskBending<T, D>::setAmbientInterface,
           py::arg("ambientInterface"), py::arg("ambientSign") = -1)
      .def(
          "setSolveBounds",
          [](OxidationMaskBending<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setSolveBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def("clearSolveBounds", &OxidationMaskBending<T, D>::clearSolveBounds)
      .def("apply", &OxidationMaskBending<T, D>::apply)
      .def("setParameters", &OxidationMaskBending<T, D>::setParameters)
      .def("getParameters", &OxidationMaskBending<T, D>::getParameters)
      .def("getIterations", &OxidationMaskBending<T, D>::getIterations)
      .def("getResidual", &OxidationMaskBending<T, D>::getResidual)
      .def("getLastApplyVelocityChange",
           &OxidationMaskBending<T, D>::getLastApplyVelocityChange)
      .def("getNumberOfSolutionNodes",
           &OxidationMaskBending<T, D>::getNumberOfSolutionNodes)
      .def("getNumberOfContactNodes",
           &OxidationMaskBending<T, D>::getNumberOfContactNodes)
      .def("getNumberOfFixedNodes",
           &OxidationMaskBending<T, D>::getNumberOfFixedNodes)
      .def("writeFieldsToLevelSet",
           &OxidationMaskBending<T, D>::writeFieldsToLevelSet);

  py::class_<OxidationConstrainedAmbient<T, D>, VelocityField<T>,
             SmartPointer<OxidationConstrainedAmbient<T, D>>>(
      module, "OxidationConstrainedAmbient")
      .def(py::init(
               [](SmartPointer<OxidationDeformation<T, D>> &deformationField,
                  SmartPointer<OxidationMaskBending<T, D>> &maskVelocityField,
                  SmartPointer<Domain<T, D>> &maskInterface, int maskSign) {
                 return OxidationConstrainedAmbient<T, D>::New(
                     deformationField, maskVelocityField, maskInterface,
                     maskSign);
               }),
           py::arg("deformationField"), py::arg("maskVelocityField"),
           py::arg("maskInterface"), py::arg("maskSign") = 1);

  py::class_<OxidationModel<T, D>, SmartPointer<OxidationModel<T, D>>>(
      module, "OxidationModel")
      .def(py::init(
               [](SmartPointer<OxidationDiffusion<T, D>> &diffusionField,
                  SmartPointer<OxidationDeformation<T, D>> &deformationField,
                  OxidationCouplingParameters<T> couplingParameters) {
                 return OxidationModel<T, D>::New(
                     diffusionField, deformationField, couplingParameters);
               }),
           py::arg("diffusionField"), py::arg("deformationField"),
           py::arg("couplingParameters") = OxidationCouplingParameters<T>())
      .def("setDiffusionField", &OxidationModel<T, D>::setDiffusionField)
      .def("setDeformationField", &OxidationModel<T, D>::setDeformationField)
      .def("setParameters", &OxidationModel<T, D>::setParameters)
      .def(
          "setSolveBounds",
          [](OxidationModel<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setSolveBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def("clearSolveBounds", &OxidationModel<T, D>::clearSolveBounds)
      .def("apply", &OxidationModel<T, D>::apply)
      .def("getIterations", &OxidationModel<T, D>::getIterations)
      .def("getResidual", &OxidationModel<T, D>::getResidual);

  py::class_<Oxidation<T, D>, SmartPointer<Oxidation<T, D>>>(module,
                                                             "Oxidation")
      .def(py::init([](SmartPointer<Domain<T, D>> &siInterface,
                       SmartPointer<Domain<T, D>> &ambientInterface,
                       SmartPointer<Domain<T, D>> &maskInterface) {
             return Oxidation<T, D>::New(siInterface, ambientInterface,
                                         maskInterface);
           }),
           py::arg("siInterface"), py::arg("ambientInterface"),
           py::arg("maskInterface"))
      .def("setSiInterface", &Oxidation<T, D>::setSiInterface)
      .def("setAmbientInterface", &Oxidation<T, D>::setAmbientInterface)
      .def("setMaskInterface", &Oxidation<T, D>::setMaskInterface)
      .def("setOxidationParameters", &Oxidation<T, D>::setOxidationParameters)
      .def("setDeformationParameters",
           &Oxidation<T, D>::setDeformationParameters)
      .def("setCouplingParameters", &Oxidation<T, D>::setCouplingParameters)
      .def("setMaskParameters", &Oxidation<T, D>::setMaskParameters)
      .def("setSpatialScheme", &Oxidation<T, D>::setSpatialScheme)
      .def("setTemporalScheme", &Oxidation<T, D>::setTemporalScheme)
      .def("setMaskCouplingIterations",
           &Oxidation<T, D>::setMaskCouplingIterations)
      .def("setMaskCouplingTolerance",
           &Oxidation<T, D>::setMaskCouplingTolerance)
      .def(
          "setSolveBounds",
          [](Oxidation<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setSolveBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def(
          "setMaskBendingBounds",
          [](Oxidation<T, D> &model,
             std::array<viennahrle::IndexType, D> passedMinIndex,
             std::array<viennahrle::IndexType, D> passedMaxIndex) {
            viennahrle::Index<D> minIndex{};
            viennahrle::Index<D> maxIndex{};
            for (unsigned i = 0; i < D; ++i) {
              minIndex[i] = passedMinIndex[i];
              maxIndex[i] = passedMaxIndex[i];
            }
            model.setMaskBendingBounds(minIndex, maxIndex);
          },
          py::arg("minIndex"), py::arg("maxIndex"))
      .def("apply", &Oxidation<T, D>::apply, py::arg("advectionTime"))
      .def("applyCFLLimited", &Oxidation<T, D>::applyCFLLimited,
           py::arg("requestedTime"), py::arg("cflFactor"),
           "Execute one CFL-limited LOCOS step; returns the actual time "
           "advanced.")
      .def("getLastMaxVelocity", &Oxidation<T, D>::getLastMaxVelocity,
           "Maximum interface velocity (µm/hr) from the most recent CFL step.")
      .def("getDiffusionField", &Oxidation<T, D>::getDiffusionField)
      .def("getDeformationField", &Oxidation<T, D>::getDeformationField)
      .def("getMaskBendingField", &Oxidation<T, D>::getMaskBendingField)
      .def("getMaskCouplingIterations",
           &Oxidation<T, D>::getMaskCouplingIterations)
      .def("getMaskCouplingResidual",
           &Oxidation<T, D>::getMaskCouplingResidual);

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
      .def(
          py::init(
              &SmartPointer<ToSurfaceMesh<T, D>>::template New<double, double>),
          py::arg("minNodeDistFactor") = 0.05, py::arg("eps") = 1e-12)
      .def(py::init(&SmartPointer<ToSurfaceMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &,
                    double, double>),
           py::arg("domain"), py::arg("mesh"),
           py::arg("minNodeDistFactor") = 0.05, py::arg("eps") = 1e-12)
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
           py::arg("minNodeDistFactor") = 0.05, py::arg("eps") = 1e-12)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    SmartPointer<Domain<T, D>> &, SmartPointer<Mesh<T>> &,
                    double, double>),
           py::arg("domain"), py::arg("mesh"),
           py::arg("minNodeDistFactor") = 0.05, py::arg("eps") = 1e-12)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    std::vector<SmartPointer<Domain<T, D>>> &,
                    SmartPointer<Mesh<T>> &, double, double>),
           py::arg("domains"), py::arg("mesh"),
           py::arg("minNodeDistFactor") = 0.05, py::arg("eps") = 1e-12)
      .def(py::init(&SmartPointer<ToMultiSurfaceMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &, double, double>),
           py::arg("mesh"), py::arg("minNodeDistFactor") = 0.05,
           py::arg("eps") = 1e-12)
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

  // ToHullMesh
  py::class_<ToHullMesh<T, D>, SmartPointer<ToHullMesh<T, D>>>(module,
                                                               "ToHullMesh")
      // constructors
      .def(py::init(&SmartPointer<ToHullMesh<T, D>>::template New<>))
      .def(py::init(&SmartPointer<ToHullMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<ToHullMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &,
                    std::vector<SmartPointer<Domain<T, D>>> &>))
      .def(py::init(&SmartPointer<ToHullMesh<T, D>>::template New<
                    SmartPointer<Mesh<T>> &, SmartPointer<Domain<T, D>> &>))
      // methods
      .def("setMesh", &ToHullMesh<T, D>::setMesh, "Set the mesh to generate.")
      .def("insertNextLevelSet", &ToHullMesh<T, D>::insertNextLevelSet,
           "Insert next level set to convert. Bigger level sets wrapping "
           "smaller ones, should be inserted last.")
      .def("clearLevelSets", &ToHullMesh<T, D>::clearLevelSets,
           "Clear all inserted level sets.")
      .def("setSharpCorners", &ToHullMesh<T, D>::setSharpCorners,
           "Set whether to generate sharp corners. Defaults to false.")
      .def("setMaterialMap", &ToHullMesh<T, D>::setMaterialMap,
           "Set the material map to use for the hull mesh.")
      .def("setBottomExtension", &ToHullMesh<T, D>::setBottomExtension,
           "Set the bottom extension value for 2D hull meshes.")
      .def("apply", &ToHullMesh<T, D>::apply, "Generate hull mesh.");

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
          [](CompareCriticalDimensions<T, D> &self, size_t index) -> py::tuple {
            T posRef{}, posCmp{}, diff{};
            bool valid =
                self.getCriticalDimensionResult(index, posRef, posCmp, diff);
            if (valid)
              return py::make_tuple(true, posRef, posCmp, diff);
            return py::make_tuple(false, T(0), T(0), T(0));
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
