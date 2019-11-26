# Contributing

### This is a check list to go through before merging any PR:

* Make sure everything builds with & without shared libs

* Run clang-format on ALL files of the project (use format-project.sh)

* Delete html folder and run make_doxygen.sh in docs/doxygen to update the website

* Wrap all implemented interface functions for Python in Wrapping/pyWrap.cpp

* IMPORTANT: Check the ReadMe file in / to make sure nothing changed
