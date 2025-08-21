/*
  This file is used to generate the python module of viennals.
  It uses pybind11 to create the modules.

  All necessary headers are included here and the interface
  of the classes which should be exposed defined
*/
#include "pyWrap.hpp"

#include <pybind11/native_enum.h>

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

  // --------- Logger ---------
  py::enum_<LogLevel>(module, "LogLevel", py::module_local())
      .value("ERROR", LogLevel::ERROR)
      .value("WARNING", LogLevel::WARNING)
      .value("INFO", LogLevel::INFO)
      .value("INTERMEDIATE", LogLevel::INTERMEDIATE)
      .value("TIMING", LogLevel::TIMING)
      .value("DEBUG", LogLevel::DEBUG);

  py::class_<Logger, SmartPointer<Logger>>(module, "Logger", py::module_local())
      .def_static("setLogLevel", &Logger::setLogLevel)
      .def_static("getLogLevel", &Logger::getLogLevel)
      .def_static("setLogFile", &Logger::setLogFile)
      .def_static("appendToLogFile", &Logger::appendToLogFile)
      .def_static("closeLogFile", &Logger::closeLogFile)
      .def_static("getInstance", &Logger::getInstance,
                  py::return_value_policy::reference)
      .def("addDebug", &Logger::addDebug)
      .def("addTiming", (Logger & (Logger::*)(const std::string &, double)) &
                            Logger::addTiming)
      .def("addTiming",
           (Logger & (Logger::*)(const std::string &, double, double)) &
               Logger::addTiming)
      .def("addInfo", &Logger::addInfo)
      .def("addWarning", &Logger::addWarning)
      .def("addError", &Logger::addError, py::arg("s"),
           py::arg("shouldAbort") = true)
      .def("print", [](Logger &instance) { instance.print(std::cout); });

  // ------ ENUMS ------
  py::native_enum<IntegrationSchemeEnum>(module, "IntegrationSchemeEnum",
                                         "enum.IntEnum")
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
             IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
      .finalize();

  py::native_enum<BooleanOperationEnum>(module, "BooleanOperationEnum",
                                        "enum.IntEnum")
      .value("INTERSECT", BooleanOperationEnum::INTERSECT)
      .value("UNION", BooleanOperationEnum::UNION)
      .value("RELATIVE_COMPLEMENT", BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .value("INVERT", BooleanOperationEnum::INVERT)
      .finalize();

  py::native_enum<CurvatureEnum>(module, "CurvatureEnum", "enum.IntEnum")
      .value("MEAN_CURVATURE", CurvatureEnum::MEAN_CURVATURE)
      .value("GAUSSIAN_CURVATURE", CurvatureEnum::GAUSSIAN_CURVATURE)
      .value("MEAN_AND_GAUSSIAN_CURVATURE",
             CurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE)
      .finalize();

  py::native_enum<FeatureDetectionEnum>(module, "FeatureDetectionEnum",
                                        "enum.IntEnum")
      .value("CURVATURE", FeatureDetectionEnum::CURVATURE)
      .value("NORMALS_ANGLE", FeatureDetectionEnum::NORMALS_ANGLE)
      .finalize();

  py::native_enum<BoundaryConditionEnum>(module, "BoundaryConditionEnum",
                                         "enum.IntEnum")
      .value("REFLECTIVE_BOUNDARY", BoundaryConditionEnum::REFLECTIVE_BOUNDARY)
      .value("INFINITE_BOUNDARY", BoundaryConditionEnum::INFINITE_BOUNDARY)
      .value("PERIODIC_BOUNDARY", BoundaryConditionEnum::PERIODIC_BOUNDARY)
      .value("POS_INFINITE_BOUNDARY",
             BoundaryConditionEnum::POS_INFINITE_BOUNDARY)
      .value("NEG_INFINITE_BOUNDARY",
             BoundaryConditionEnum::NEG_INFINITE_BOUNDARY)
      .finalize();

  py::native_enum<FileFormatEnum>(module, "FileFormatEnum", "enum.IntEnum")
      .value("VTK_LEGACY", FileFormatEnum::VTK_LEGACY)
      .value("VTP", FileFormatEnum::VTP)
      .value("VTU", FileFormatEnum::VTU)
      .finalize();

  py::native_enum<VoidTopSurfaceEnum>(module, "VoidTopSurfaceEnum",
                                      "enum.IntEnum")
      .value("LEX_LOWEST", VoidTopSurfaceEnum::LEX_LOWEST)
      .value("LEX_HIGHEST", VoidTopSurfaceEnum::LEX_HIGHEST)
      .value("LARGEST", VoidTopSurfaceEnum::LARGEST)
      .value("SMALLEST", VoidTopSurfaceEnum::SMALLEST)
      .finalize();

  py::native_enum<TransformEnum>(module, "TransformEnum", "enum.IntEnum")
      .value("TRANSLATION", TransformEnum::TRANSLATION)
      .value("ROTATION", TransformEnum::ROTATION)
      .value("SCALE", TransformEnum::SCALE)
      .finalize();


  // ---------- MESH CLASSES ----------
  // Mesh
  py::class_<Mesh<T>, SmartPointer<Mesh<T>>>(module, "Mesh")
      // constructors
      .def(py::init(&SmartPointer<Mesh<T>>::New<>))
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

  // TransformMesh
  py::class_<TransformMesh<T>, SmartPointer<TransformMesh<T>>>(module,
                                                               "TransformMesh")
      // constructors
      .def(py::init([](SmartPointer<Mesh<T>> &mesh, TransformEnum op,
                       Vec3D<T> vec, double angle) {
             return SmartPointer<TransformMesh<T>>::New(mesh, op, vec, angle);
           }),
           py::arg("mesh"), py::arg("transform") = TransformEnum::TRANSLATION,
           py::arg("transformVector") = Vec3D<T>{0., 0., 0.},
           py::arg("angle") = 0.)
      // methods
      .def("apply", &TransformMesh<T>::apply, "Apply the transformation.");

  // VTKReader
  py::class_<VTKReader<T>, SmartPointer<VTKReader<T>>>(module, "VTKReader")
      // constructors
      .def(py::init(&SmartPointer<VTKReader<T>>::New<>))
      .def(py::init(&SmartPointer<VTKReader<T>>::New<SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<VTKReader<T>>::New<SmartPointer<Mesh<T>> &,
                                                     std::string>))
      .def(py::init([](SmartPointer<Mesh<T>> &mesh, FileFormatEnum format,
                       std::string s) {
        return SmartPointer<VTKReader<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &VTKReader<T>::setMesh, "Set the mesh to read into.")
      .def("setFileFormat", &VTKReader<T>::setFileFormat,
           "Set the file format of the file to be read.")
      .def("setFileName", &VTKReader<T>::setFileName,
           "Set the name of the input file.")
      .def("getMetaData", &VTKReader<T>::getMetaData,
           "Get the metadata from the file.")
      .def("apply", &VTKReader<T>::apply, "Read the mesh.");

  // VTKWriter
  py::class_<VTKWriter<T>, SmartPointer<VTKWriter<T>>>(module, "VTKWriter")
      // constructors
      .def(py::init(&SmartPointer<VTKWriter<T>>::New<>))
      .def(py::init(&SmartPointer<VTKWriter<T>>::New<SmartPointer<Mesh<T>> &>))
      .def(py::init(&SmartPointer<VTKWriter<T>>::New<SmartPointer<Mesh<T>> &,
                                                     std::string>))
      .def(py::init([](SmartPointer<Mesh<T>> &mesh, FileFormatEnum format,
                       std::string s) {
        return SmartPointer<VTKWriter<T>>::New(mesh, format, s);
      }))
      // methods
      .def("setMesh", &VTKWriter<T>::setMesh, "Set the mesh to output.")
      .def("setFileFormat", &VTKWriter<T>::setFileFormat,
           "Set the file format, the mesh should be written to.")
      .def("setFileName", &VTKWriter<T>::setFileName,
           "Set the name of the output file.")
      .def("setMetaData", &VTKWriter<T>::setMetaData,
           "Set the metadata to be written to the file.")
      .def(
          "addMetaData",
          py::overload_cast<const std::string &, T>(&VTKWriter<T>::addMetaData),
          "Add a single metadata entry to the file.")
      .def("addMetaData",
           py::overload_cast<const std::string &, const std::vector<T> &>(
               &VTKWriter<T>::addMetaData),
           "Add a single metadata entry to the file.")
      .def("addMetaData",
           py::overload_cast<
               const std::unordered_map<std::string, std::vector<T>> &>(
               &VTKWriter<T>::addMetaData),
           "Add metadata to the file.")
      .def("apply", &VTKWriter<T>::apply, "Write the mesh.");

  // --------------- OTHER CLASSES ----------------
  // MaterialMap
  py::class_<MaterialMap, SmartPointer<MaterialMap>>(module, "MaterialMap")
      // constructors
      .def(py::init(&SmartPointer<MaterialMap>::New<>))
      // methods
      .def("insertNextMaterial", &MaterialMap::insertNextMaterial,
           "Insert a new material into the map.")
      .def("setMaterialId", &MaterialMap::setMaterialId)
      .def("getNumberOfLayers", &MaterialMap::getNumberOfLayers,
           "Get the number of level-sets in the material map.")
      .def("getNumberOfMaterials", &MaterialMap::getNumberOfMaterials)
      .def("getMaterialId", &MaterialMap::getMaterialId);

  // PointData
  py::class_<PointData<T>, SmartPointer<PointData<T>>>(module, "PointData")
      // constructors
      .def(py::init(&SmartPointer<PointData<T>>::New<>))
      // methods
      .def("insertNextScalarData",
           (void(PointData<T>::*)(const PointData<T>::ScalarDataType &,
                                  const std::string &)) &
               PointData<T>::insertNextScalarData,
           py::arg("scalars"), py::arg("label") = "Scalars")
      .def("insertNextVectorData",
           (void(PointData<T>::*)(const PointData<T>::VectorDataType &,
                                  const std::string &)) &
               PointData<T>::insertNextVectorData,
           py::arg("vectors"), py::arg("label") = "Vectors")
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

  // Slice
  py::class_<Slice<T>, SmartPointer<Slice<T>>>(module, "Slice")
      // constructors
      .def(py::init(&SmartPointer<Slice<T>>::New<>))
      .def(py::init(
          &SmartPointer<Slice<T>>::New<SmartPointer<Domain<T, 3>> &,
                                       SmartPointer<Domain<T, 2>> &, int, T>))
      .def(py::init(
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

  // Extrude
  py::class_<Extrude<T>, SmartPointer<Extrude<T>>>(module, "Extrude")
      // constructors
      .def(py::init(&SmartPointer<Extrude<T>>::New<>))
      .def(py::init(
          &SmartPointer<Extrude<T>>::New<SmartPointer<Domain<T, 2>> &,
                                         SmartPointer<Domain<T, 3>> &,
                                         std::array<T, 2>, const int,
                                         std::array<BoundaryConditionEnum, 3>>))
      // methods
      .def("setInputLevelSet", &Extrude<T>::setInputLevelSet,
           "Set 2D input Level Set")
      .def("setOutputLevelSet", &Extrude<T>::setOutputLevelSet,
           "Set 3D output Level Set")
      .def("setExtent", &Extrude<T>::setExtent,
           "Set the extent in the extruded dimension")
      .def("setExtrudeDimension", &Extrude<T>::setExtrudeDimension,
           "Set the dimension which should be extruded")
      .def("setBoundaryConditions",
           py::overload_cast<std::array<BoundaryConditionEnum, 3>>(
               &Extrude<T>::setBoundaryConditions),
           "Set the boundary conditions in the 3D extruded domain.")
      .def("setBoundaryConditions",
           py::overload_cast<BoundaryConditionEnum *>(
               &Extrude<T>::setBoundaryConditions),
           "Set the boundary conditions in the 3D extruded domain.")
      .def("apply", &Extrude<T>::apply, "Perform extrusion.");

  // VelocityField
  py::class_<VelocityField<T>, SmartPointer<VelocityField<T>>,
             PylsVelocityField>(module, "VelocityField")
      // constructors
      .def(py::init<>())
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

  // ---------- MAIN API ----------
  // Submodule for 2D
  auto m2 = module.def_submodule("d2", "2D bindings");
  bindApi<2>(m2);

  // Submodule for 3D
  auto m3 = module.def_submodule("d3", "3D bindings");
  bindApi<3>(m3);
}
