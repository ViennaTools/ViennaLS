# ViennaLS

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

```
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
make buildDependencies # this will install all dependencies and might take a while
make install
```

This will install the necessary headers and CMake files to the specified path. If DCMAKE_INSTALL_PREFIX is not specified, it will be installed to the standard path for your system, usually /usr/local/ .

## Installing without VTK

In order to install ViennaLS without VTK, run:
```
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -DVIENNALS_USE_VTK=OFF
make buildDependencies
make install
```

## Installing with dependencies already installed on the system

If you want to use your own install of dependencies, just specify the directories of dependencies in CMake:

```bash
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -DVTK_DIR=/path/to/vtk/install -Dpybind11_DIR=/path/to/pybind11 -DViennaHRLE_DIR=/path/to/viennahrle
make install
```

It is possible to specify one preinstalled dependency and automatically install the others, by just passing the path to one dependency as shown above, but not the others, i.e. for using preinstalled vtk but automatically installed ViennaHRLE and pybind11:

```bash
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/ -DVTK_DIR=/path/to/vtk/install
make buildDependencies
make install
```

## Using the viennaLS python module

The Releases only contain the compiled library for the most common Python version per platform:
* Windows: Python 3.8
* Linux: Python 3.8

For all other Python versions, you have to build the library yourself (see below).

In order to use ViennaLS in python, just download the python shared libraries from the [releases section](https://github.com/ViennaTools/ViennaLS/releases) and put it in your current folder.
From this folder just import the 2D or the 3D version of the library:

```
import viennaLS2d as vls
levelset = vls.lsDomain(0.2) # empty level set with grid spacing 0.2
sphere = vls.lsSphere((0,0,0), 5) # sphere at origin with radius 5
vls.lsMakeGeometry(levelset, sphere).apply() # create sphere in level set
```

All functions which are available in C++ are also available in Python. In order to switch to three dimensions, only the import needs to be changed:

```
import viennaLS3d as vls
```


## Building the python module

In order to build the python module, set `VIENNALS_BUILD_PYTHON_2` or `VIENNALS_BUILD_PYTHON_3` to `ON`:
```
cmake .. -DVIENNALS_BUILD_PYTHON_3=ON
make buildDependencies # this will install pybind11 the first time it is called
make
```

If both options are on, only VIENNALS_BUILD_PYTHON_3 will be used, since only one version can be built at a time.

## Setting up the dependencies

If you just want to install all dependencies before doing anything else, run:

```
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
make pybind11-external
make viennahrle-external
make vtk-external
```

This will take some time the first time it is run.
The dependencies only need to be set up once.

## Running the Tests

ViennaLS uses CTest to run its tests.
In order to check whether ViennaLS runs without issues on your system, you can run:

```
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DVIENNALS_BUILD_TESTS=ON
make buildTests # build all tests
make test # run all tests
```

## Building examples

The examples can be built using CMake:

```
cmake .. -DVIENNALS_BUILD_EXAMPLES=ON
make
```

## Integration in CMake projects

In order to use this library in your CMake project, add the following lines to the CMakeLists.txt of your project:

```
set(ViennaLS_DIR "/path/to/your/custom/install/")
find_package(ViennaLS REQUIRED PATHS ${ViennaLS_DIR})
add_executable(myExe mySource.cpp)
target_include_directories(myExe PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(myExe ${VIENNALS_LIBRARIES})
```

### Shared libraries

In order to save build time during developement, dynamically linked shared libraries can be used
if ViennaLS was built with them. This is done by precompiling the most common template specialisations.
In order to use shared libraries, use 
```
cmake .. -DVIENNALS_BUILD_SHARED_LIBS=ON
```
If ViennaLS was build with shared libraries and you use ViennaLS in your project (see above), CMake will automatically link them to your project. In order to build a release of your own project with better runtime performance, but
longer build times, use the following CMake option when building a release:
```
VIENNALS_USE_SHARED_LIBS=OFF
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
