<div align="center">

![logo](assets/logo.png)

# ViennaLS

[![Linux](https://github.com/ViennaTools/ViennaLS/actions/workflows/linux_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/linux_test.yml)
[![macOS](https://github.com/ViennaTools/ViennaLS/actions/workflows/macos_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/macos_test.yml)
[![Windows](https://github.com/ViennaTools/ViennaLS/actions/workflows/windows_test.yml/badge.svg)](https://github.com/ViennaTools/ViennaLS/actions/workflows/windows_test.yml)

</div>

ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualisation applications, although this is not the main design target.

> [!NOTE]  
> ViennaLS is under heavy development and improved daily. If you do have suggestions or find bugs, please let us know!

## Support

[Documentation](https://viennatools.github.io/ViennaLS/index.html) and [Examples](https://github.com/ViennaTools/ViennaLS/tree/master/Examples) can be found online.

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

The CMake configuration automatically checks if the dependencies are installed. If CMake is unable to find them, the dependencies will be built from source.

## Building the Python package

> [!NOTE]  
> On systems that feature a package manager (e.g. Ubuntu/Debian `apt`), the dependencies can be installed beforehand (e.g. using ```sudo apt install libvtk9.1 libvtk9-dev pybind11-dev```), which saves a considerable amount of time during compilation.

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
  CPMAddPackage("gh:viennatools/viennals@3.0.0")
  ```

* With a local installation
    > In case you have ViennaLS installed in a custom directory, make sure to properly specify the `CMAKE_MODULE_PATH`
    ```cmake
    find_package(ViennaLS REQUIRED)
    target_link_libraries(${PROJECT_NAME} ViennaTools::ViennaLS)
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

If you want to contribute to ViennaLS, make sure to follow the [LLVM Coding guidelines](https://llvm.org/docs/CodingStandards.html). Before creating a pull request, make sure ALL files have been formatted by clang-format, which can be done using the `format-project.sh` script in the root directory.

## Authors

Current contributors: Lado Filipovic, Paul Manstetten, Xaver Klemenschits and Josef Weinbub

Founder and initial developer: Otmar Ertl

Contact us via: viennats@iue.tuwien.ac.at

ViennaLS was developed under the aegis of the 'Institute for Microelectronics' at the 'TU Wien'.
http://www.iue.tuwien.ac.at/

## License

See file LICENSE in the base directory.
