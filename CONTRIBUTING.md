# Contributing

### This is a check list to go through before merging any PR:

* Make sure everything builds with & without shared libs

* Run clang-format on ALL files of the project (use CMake target `format`)

* Wrap all implemented interface functions for Python in python/pyWrap.cpp

* IMPORTANT: Check the ReadMe file in / to make sure nothing changed
