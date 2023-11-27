# ViennaLS
[![Linux](https://github.com/ViennaTools/ViennaLS/actions/workflows/linux_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/linux_test.yml)
[![macOS](https://github.com/ViennaTools/ViennaLS/actions/workflows/macos_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/macos_test.yml)
[![Windows](https://github.com/ViennaTools/ViennaLS/actions/workflows/windows_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/windows_test.yml)

ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualisation applications, although this is not the main design target.

IMPORTANT NOTE: ViennaLS is under heavy development and improved daily. If you do have suggestions or find bugs, please let us know!

## Support

[Documentation](https://viennatools.github.io/ViennaLS/doxygen/html/index.html) and [Examples](https://viennatools.github.io/ViennaLS/doxygen/html/examples.html) can be found online.

Bug reports and suggestions should be filed on GitHub.

## Releases
Releases are tagged on the maser branch and available in the [releases section](https://github.com/ViennaTools/ViennaLS/releases).

## Building

### Supported Operating Systems

* Windows (Visual Studio)

* Linux (g++ / clang)

* macOS (XCode)


### System Requirements

* C++17 Compiler with OpenMP support

### Dependencies (installed automatically)

* [ViennaHRLE](https://github.com/ViennaTools/ViennaHRLE)

* [VTK](https://github.com/Kitware/VTK) (optional)

* [pybind11](https://github.com/pybind/pybind11) (only for building Python libs)

## Using ViennaLS in your project

Have a look at the [example repo](https://github.com/ViennaTools/viennals-example) for creating a project with ViennaLS as a dependency.


## Installing

Since this is a header only project, it does not require any installation.
However, we recommend the following procedure in order to set up all dependencies correctly:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --build build
```

This will install the necessary headers and CMake files to the specified path. If DCMAKE_INSTALL_PREFIX is not specified, it will be installed to the standard path for your system, usually /usr/local/ .

## Installing without VTK

In order to install ViennaLS without VTK, run:
```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -DVIENNALS_USE_VTK=OFF
cmake --build build
```

## Installing with dependencies already installed on the system

The CMake configuration automatically checks if the dependencies are installed. If CMake is unable to find them, the dependencies will be built from source with the _buildDependencies_ target.

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
cmake --build build
```
If one wants to use a specific installation of one or more of the dependencies, just pass the corresponding _*_DIR_ variable as a configuration option (e.g. -DVTK_DIR=/path/to/vtk/install -Dpybind11_DIR=/path/to/pybind11 -DViennaHRLE_DIR=/path/to/viennahrle)

## Building the Python package

> __Tip__: On systems that feature a package manager (e.g. Ubuntu/Debian `apt`), the dependencies can be installed beforehand (e.g. using ```sudo apt install libvtk9.1 libvtk9-dev pybind11-dev```), which saves a considerable amount of time during compilation.

The Python package can be built and installed using the `pip` command:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS
pip install --user .
```

## Using the Python package

All functions which are available in C++ are also available in Python. The 2D version of the library can be imported as follows:
```python
import viennals2d as vls
```

In order to switch to three dimensions, only the import needs to be changed:

```python
import viennals3d as vls
```

## Setting up the dependencies

If you just want to install all dependencies before doing anything else, run:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
make pybind11_external
make viennahrle_external
make vtk_external
```

This will take some time the first time it is run.
The dependencies only need to be set up once.

## Running the Tests

ViennaLS uses CTest to run its tests.
In order to check whether ViennaLS runs without issues on your system, you can run:

```bash
git clone https://github.com/ViennaTools/ViennaLS.git
cd ViennaLS

cmake -B build -DVIENNALS_BUILD_TESTS=ON
cmake --build build
ctest -LE '^benchmark$' --test-dir build
```

## Building examples

The examples can be built using CMake:

```bash
cmake -B build -DVIENNALS_BUILD_EXAMPLES=ON
cmake --build build
```

## Integration in CMake projects

In order to use this library in your CMake project, add the following lines to the CMakeLists.txt of your project:

```cmake
set(ViennaLS_DIR "/path/to/your/custom/install/")
find_package(ViennaLS REQUIRED PATHS ${ViennaLS_DIR})
add_executable(myExe mySource.cpp)
target_include_directories(myExe PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(myExe ${VIENNALS_LIBRARIES})
```

### Shared libraries

In order to save build time during development, dynamically linked shared libraries can be used
if ViennaLS was built with them. This is done by precompiling the most common template specialisations.
In order to use shared libraries, use 
```bash
cmake .. -DVIENNALS_PRECOMPILE_HEADERS=ON
```
If ViennaLS was built with shared libraries and you use ViennaLS in your project (see above), CMake will automatically link them to your project. In order to build a release of your own project with better runtime performance, but
longer build times, use the following CMake option when building a release:
```bash
VIENNALS_USE_PRECOMPILED=OFF
```

## Contributing
Before being able to merge your PR, make sure you have met all points on the checklist in [CONTRIBUTING.md](https://github.com/ViennaTools/viennals/blob/master/CONTRIBUTING.md).

If you want to contribute to ViennaLS, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html). Before creating a pull request, make sure ALL files have been formatted by clang-format, which can be done using the format-project.sh script in the root directory.

## Authors

Current contributors: Lado Filipovic, Paul Manstetten, Xaver Klemenschits and Josef Weinbub

Founder and initial developer: Otmar Ertl

Contact us via: viennats@iue.tuwien.ac.at

ViennaLS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

License
--------------------------
See file LICENSE in the base directory.
