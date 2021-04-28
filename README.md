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

* [ViennaHRLE](https://github.com/ViennaTools/viennahrle)

* [VTK](https://vtk.org/) (optional)


## Using ViennaLS in your project

Have a look at the [example repo](https://github.com/ViennaTools/viennals-example) for creating a project with ViennaLS as a dependency.


## Installing (with dependencies already installed)

Since this is a header only project, it does not require any installation.
However, we recommend the following procedure.

Make sure you have [ViennaHRLE](https://github.com/ViennaTools/viennahrle) installed on your system and run:

```
git clone github.com/ViennaTools/ViennaLS.git
cd ViennaLS
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/custom/install/
make install
```

This will install the necessary headers and CMake files to the specified path. If DCMAKE_INSTALL_PREFIX is not specified, it will be installed to the standard path for your system, usually /usr/local/ .


## Using the viennaLS python module

The Releases only contain the compiled library for the most common Python version per platform:
* Windows: Python 3.8
* Linux: Python 3.6

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


### Building the python module

In order to build the python module, set `VIENNALS_BUILD_PYTHON_2_7` or `VIENNALS_BUILD_PYTHON_3_6` to `ON`:
```
cmake .. -DVIENNALS_BUILD_PYTHON_3_6=ON
make
```

If both options are on, only VIENNALS_BUILD_PYTHON_3_6 will be used, since only one version can be built at a time.

## Integration in CMake projects

In order to use this library in your CMake project, add the following lines to the CMakeLists.txt of your project:\
(also do not forget to include ViennaHRLE)

```
set(ViennaLS_DIR "/path/to/your/custom/install/")
find_package(ViennaLS REQUIRED)
add_executable(...)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${VIENNALS_LIBRARIES})
```

### Building examples

The examples can be built using CMake:

```
mkdir build && cd build
cmake .. -DVIENNALS_BUILD_EXAMPLES=ON
make
```

### Shared libraries

In order to save build time during developement, dynamically linked shared libraries can be used
if ViennaLS was built with them. This is done by precompiling the most common template specialisations.
In order to use shared libraries, use 
```
cmake .. -DVIENNALS_BUILD_SHARED_LIBS=ON
```
If ViennaLS was build with shared libraries and you use ViennaLS in your project (see above), CMake will automatically link them to your project. In order to build a release of your own project with better runtime performance, but
longer build times, use the following CMake option when building your project:
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
