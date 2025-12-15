<div align="center">

<img src="assets/logo.svg" height="100" />
<h1>ViennaLS</h1>

[![🧪 Tests](https://github.com/ViennaTools/ViennaLS/actions/workflows/build.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/build.yml)
[![🐍 Bindings](https://github.com/ViennaTools/ViennaLS/actions/workflows/python.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/python.yml)
[![PyPi Version](https://img.shields.io/pypi/v/ViennaLS?logo=pypi)](https://pypi.org/project/ViennaLS/)

</div>

ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualisation applications, although this is not the main design target.

> [!NOTE]  
> ViennaLS is under heavy development and improved daily. If you do have suggestions or find bugs, please let us know!

## Quick Start  

To install ViennaLS for Python, simply run:  

```sh
pip install ViennaLS
```

To use ViennaLS in C++, clone the repository and follow the installation steps below.

## Support

[Documentation](https://viennatools.github.io/ViennaLS/index.html) and [Examples](https://github.com/ViennaTools/ViennaLS/tree/master/examples) can be found online.

Bug reports and suggestions should be filed on GitHub.

## Releases

Releases are tagged on the maser branch and available in the [releases section](https://github.com/ViennaTools/ViennaLS/releases).

## Building

### Supported Operating Systems

* Windows (MSVC)

* Linux (g++ & clang)

* macOS (XCode)


### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies

> Dependencies will be installed automatically when not available.

* [ViennaHRLE](https://github.com/ViennaTools/ViennaHRLE)

* [VTK](https://github.com/Kitware/VTK) (optional, but recommended for mesh export and visualization)

* [pybind11](https://github.com/pybind/pybind11) (only for building Python libs)

## Using ViennaLS in your project

Have a look at the [example repo](https://github.com/ViennaTools/viennals-example) for creating a project with ViennaLS as a dependency.

## Installing

Since this is a header only project, it does not require any installation.
However, we recommend the following procedure in order to set up all dependencies correctly:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -D CMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --install build
```

This will install the necessary headers and CMake files to the specified path. If `CMAKE_INSTALL_PREFIX` is not specified, it will be installed to the standard path for your system, usually `/usr/local/`.

## Installing without VTK

In order to install ViennaLS without VTK, run:
```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -D CMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -D VIENNALS_USE_VTK=OFF
cmake --install build
```

## Installing with dependencies already installed on the system

The CMake configuration automatically checks if the dependencies are installed. If CMake is unable to find them, the dependencies will be built from source.

## Building the Python package

> [!NOTE]  
> On systems that feature a package manager (e.g. Ubuntu/Debian `apt`), VTK can be installed beforehand (e.g. using ```sudo apt install libvtk9-dev```), which saves a considerable amount of time during compilation.

The Python package can be built and installed using the `pip` command:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

pip install .
```

## Using the Python package

All functions which are available in C++ are also available in Python. The 2D version of the library can be imported as follows:
```python
import viennals.d2 as vls 
import viennals # for common functions
```

To switch to 3D, only the import changes:

```python
import viennals.d3 as vls
```
Functions that operate on a domain object (e.g. `Advect`, `ToSurfaceMesh`, ...) are provided in the respective `d2` or `d3` modules.
Common functions, enums, and dimension-independent utilities (such as `Mesh`) are available directly in the `viennals` namespace.

A complete list of functions and their locations can be found in the [API documentation](PythonAPI.md).

For examples on how to use the Python package, please have a look at these examples: [Air Gap Deposition](https://github.com/ViennaTools/ViennaLS/blob/master/examples/AirGapDeposition/AirGapDeposition.py), [Deposition](https://github.com/ViennaTools/ViennaLS/blob/master/examples/Deposition/Deposition.py), [Geometric Advection](https://github.com/ViennaTools/ViennaLS/blob/master/examples/GeometricAdvection/GeometricAdvection.py).


## Running the Tests

ViennaLS uses CTest to run its tests.
In order to check whether ViennaLS runs without issues on your system, you can run:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -DVIENNALS_BUILD_TESTS=ON
cmake --build build
ctest -E "Benchmark|Performance" --test-dir build
```

## Building examples

The examples can be built using CMake:

```bash
cmake -B build -DVIENNALS_BUILD_EXAMPLES=ON
cmake --build build
```

## Integration in CMake projects

We recommend using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) to consume this library.

* Installation with CPM
  ```cmake
  CPMAddPackage("gh:viennatools/viennals@5.3.0")
  ```

* With a local installation
    > In case you have ViennaLS installed in a custom directory, make sure to properly specify the `CMAKE_MODULE_PATH` or `PATHS` in your `find_package` call.

    ```cmake
    set(VIENNALS_PATH "/your/local/installation")

    find_package(OpenMP REQUIRED)
    find_package(VTK        PATHS ${VIENNALS_PATH})
    find_package(ViennaHRLE PATHS ${VIENNALS_PATH})
    find_package(ViennaLS   PATHS ${VIENNALS_PATH})

    target_link_libraries(${PROJECT_NAME} PUBLIC ViennaTools::ViennaLS)
    ```

### Shared Library

In order to save build time during development, dynamically linked shared libraries can be used
if ViennaLS was built with them. This is done by precompiling the most common template specialisations.
In order to use shared libraries, use 
```bash
cmake -B build -DVIENNALS_PRECOMPILE_HEADERS=ON
```
If ViennaLS was built with shared libraries and you use ViennaLS in your project (see above), CMake will automatically link them to your project.

## Contributing

Before being able to merge your PR, make sure you have met all points on the checklist in [CONTRIBUTING.md](https://github.com/ViennaTools/viennals/blob/master/CONTRIBUTING.md).

If you want to contribute to ViennaLS, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html).

Make sure to format all files before creating a pull request:
```bash
cmake -B build
cmake --build build --target format
```

## Authors

Current contributors: Tobias Reiter, Roman Kostal, Lado Filipovic

Founder and initial developer: Otmar Ertl

Contact us via: viennatools@iue.tuwien.ac.at

ViennaLS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.  
http://www.iue.tuwien.ac.at/

## License

ViennaLS is licensed under the [MIT License](./LICENSE).

Some third-party libraries used by ViennaLS are under their own permissive licenses (MIT, BSD).  
See [`THIRD_PARTY_LICENSES.md`](./THIRD_PARTY_LICENSES.md) for details.
