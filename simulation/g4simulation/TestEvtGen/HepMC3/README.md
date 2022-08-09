# Introduction

HepMC3 is a new version of the HepMC event record. It uses shared pointers for in-memory navigation and the POD concept for persistency.
Visit the home page of the project http://hepmc.web.cern.ch/hepmc/ or the CERN GitLab repository https://gitlab.cern.ch/hepmc/HepMC3 for more information.
You can also send bug reports, feature requests and questions about HepMC3 to hepmc-devATcern.ch


# Quick-start (applicable to recent versions):

0. On the Linux platform HepMC3 is available from the standard repositories of multiple many
linux distributions:  Fedora, CentOS(EPEL),  Arch Linux, Mageia, openSUSE Tumbleweed, Debian (outdated) etc.
It is recommended to use the system package manager to install HepMC3 for these systems.
The following commands will install the HepMC3 in the corresponding systems:
  - Fedora:
    The package is available from the standard repository. To install:

    ```
    sudo dnf install HepMC3 HepMC3-devel HepMC3-search HepMC3-search-devel HepMC3-interfaces-devel HepMC3-doc

    ```
    To have a full installation of the HepMC3-doc package, add an option `--setopt=tsflags=''`

  - RHEL, CentOS and compatible
    The package is available from the EPEL repository. To install:

    ```
    sudo yum install epel-release
    sudo yum install HepMC3 HepMC3-devel HepMC3-search HepMC3-search-devel HepMC3-interfaces-devel HepMC3-doc

    ```
  For the Fedora and CentOS the ROOT interface packages ``HepMC3-rootIO`` and ``HepMC3-rootIO-devel``
  can be installed in the same way, but will bring root as a dependency.
  The python binding packages could have different names depending on the used system, e.g. these are
  ``python2-HepMC3``, ``python34-HepMC3``, ``python36-HepMC3``, ``python2-HepMC3-rootIO``, ``python34-HepMC3-rootIO``
  and ``python36-HepMC3-rootIO`` for the CentOS7. Therefore, the simplest solution is to install bindings for all the
  available python version in the system using
  ```
    sudo yum install python*-HepMC3
  ```
  - openSUSE/Leap
    The HepMC3 package is available from the standard repositories https://build.opensuse.org/package/show/openSUSE:Leap:15.2:Update/HepMC.
    To install:
    ```
    sudo zypper install HepMC3
    ```
    The package does not include ROOT interface.

  - Arch linux
    The package is available from https://www.archlinux.org/packages/community/x86_64/hepmc/. To install:
    ```
    sudo pacman -Syu hepmc
    ```
    The dependencies can vary.

  - MacOSX
    The HepMC3 package is available in the ``homebrew-hep`` repository https://davidchall.github.io/homebrew-hep/.
    To install:
    ```
    brew tap davidchall/hep
    brew install hepmc3
    ```
    The package optionally includes ROOT interface.

   - Windows
    The precompilled HepMC3 packages are available in for Windows and other platforms also in PyPi.
    The Windows users can use ``pip`` to install HepMC3.
    ```
    pip install HepMC3
    ```
    The packages from pip do not include ROOT interface.

   - LCG platforms
     The HepMC3 package is available in the LCG from multiple platforms.
     The HepMC3 cmake files are located under
     ``/cvmfs/sft.cern.ch/lcg/views/LCG_97/${LCG_PLATFORM}/share/HepMC3/cmake``, e.g. in
     ``/cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/share/HepMC3/cmake`` for the ``x86_64-centos7-gcc8-opt`` platform.

If for any reason the instalation from the repositories is not possible or is not desired, the instructions to build
HepMC3 from the sources
The minimal requirements for the  installation from sources are
 - C++ compiler with c++11 standard support.
 - CMake version 3.X.

It is strongly recommended to read this documentation completely before the installation.
However, if for some reason that is not possible, here is a  set of commands for the
instalation that can be copied and pasted into a Unix terminal. In some cases this
action might result in a fully functional installation.

```
  wget http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.2.6.tar.gz
  tar -xzf HepMC3-3.2.6.tar.gz
  mkdir hepmc3-build
  cd hepmc3-build
  cmake -DCMAKE_INSTALL_PREFIX=../hepmc3-install   \
        -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF            \
        -DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF        \
        -DHEPMC3_ENABLE_TEST:BOOL=OFF              \
        -DHEPMC3_INSTALL_INTERFACES:BOOL=ON        \
        -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF        \
        -DHEPMC3_BUILD_DOCS:BOOL=OFF     \
        -DHEPMC3_ENABLE_PYTHON:BOOL=ON   \
        -DHEPMC3_PYTHON_VERSIONS=2.7     \
        -DHEPMC3_Python_SITEARCH27=../hepmc3-install/lib/python2.7/site-packages \
        ../HepMC3-3.2.6
  make
  make install
```
The full explanation  for the installation options is given below.


1. The first step of the installation is to
  a)Checkout the HepMC from GIT repository:

  ```
  git clone https://gitlab.cern.ch/hepmc/HepMC3.git
  ```
  or b) download the tarball with wget or any other tool
  ```
  wget http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.2.6.tar.gz
  tar -xzf HepMC3-3.2.6.tar.gz
  ```


2. Create a workspace area on which to perform the builds
  ```
  mkdir hepmc3-build
  cd hepmc3-build
  ```


3. Configure the HepMC3 features to build and install:
  a)
  ```
  cmake -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=../hepmc3-install ../HepMC3
  ```
  or b)
  ```
  cmake -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=../hepmc3-install ../HepMC3-3.2.6
  ```

  To disable compilation of the search engine (enabled by default) add
  ```
  -DHEPMC3_ENABLE_SEARCH=OFF
  ```
  Please note, in case of using the search engine, the easiest way to link against
  HepMC3search could be to use static version of the library.

  To disable compilation of Python bindings (enabled by default) add
  ```
  -DHEPMC3_ENABLE_PYTHON=OFF
  ```
  See points 10a and 10b below for more details on configuring the Python extension
  for installation into non-system areas (e.g. in a virtual environment) and for
  compatibility with PyPy.

  Please note that on some RedHat-like systems the executable of CMake version 3 is called `cmake3`.


4. In order to build with ROOTIO put the following flags and define LD_LIBRARY_PATH:
  ```
  -DHEPMC3_ENABLE_ROOTIO=ON -DROOT_DIR=path_to_ROOT6_installation
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:path_to_ROOT6_libraries
  ```


5. In order to build with HepMC3 example put
  ```
  -DHEPMC3_BUILD_EXAMPLES=ON
  ```

  Additional flags control the location installation of Photos, Tauola MC-TESTER, HepMC2
  and Pythia8  that are used for tests. Please note that some examples are disabled  on Windows.

  For tests and for development purposes  HepMC3 uses extra packages, e.g. GrpahViz, binder.
  On Windows, most extra packages have no default location.


6. In order to build and run tests put
  ```
  -DHEPMC3_ENABLE_TEST=ON
  ```
  and run
  ```
  make test
  ```
  after compilation is completed


7. Build and install HepMC3
  ```
  make -jN install

  ```
where N denotes a number of desired parallel jobs.
On Windows (in Unix-compatible shell) one can use
```
cmake --build ./
```
instead

Please note that you have to have permissions to install HepMC3 into the selected directory.


8. To run the tests do
```
ctest .
```
On Windows/MSVC the build type should be specified
```
ctest . -C Debug
```


9. HepMC3 packages are available from standard repositories of Debian,
Redhat/Fedora, Arch, hep-homebrew and others.
To build Windows installer use NSIS and run
```
cpack.exe -G NSIS -C Debug
```
The primary supported development chains are
gcc@Linux, clang@Linux, gcc@OSX, clang@OSX and  MSVC@Windows10, however
the code has been tested on other compilers and systems, e.g. clang@FreeBSD,
Sun@Solaris11, gcc@Solaris11, Sun@Linux, PGC@Linux etc.
Only a limiter support can be offered for these systems.



10. HepMC3 is supplied with Python bindings codes that can be build on multiple systems.
The number of potential combinations of compiler suites, Python versions and operation systems is quite large,
therefore it is recommended to install the HepMC3 form the default repositories of the used operating system.
If not possible, the following options are available
- compile the Python package from sources
- use the precompiled binaries from pypi: https://pypi.org/project/HepMC3/
- use the https://anaconda.org/conda-forge/hepmc3
In addition to that part of the bindings are implemented in the pyhepmc project

https://github.com/HDembinski/pyhepmc.

Please note that


10a

The installation path for the Python modules can be tweaked with
a set of dynamically named CMake variables:
```
HEPMC3_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR}
```
For example, the variable is called `HEPMC3_Python_SITEARCH36` for Python version 3.6.
So, to set the installation path of the HepMC3 modules for Python 3.6, add this to the
arguments of `cmake`:
```
-DHEPMC3_Python_SITEARCH36=/path/where/I/want/my/python/of/version/3.6/modules
```

By default the `HEPMC3_Python_SITEARCH${Python_VERSION_MAJOR}${Python_VERSION_MINOR}`
variables are set to the path of the used installation of Python. In some cases you
would not have permissions to install the Python modules to that directory, so the change
of the installation path for the Python modules would be mandatory. A natural choice is
the Python module directory inside your installation prefix, e.g.
```
cmake ...
   -DCMAKE_INSTALL_PREFIX=/my/install/path
   -DHEPMC3_Python_SITEARCH36=/my/install/path/lib/python3.6/site-packages
```

To specify the desired versions of Python to build the HepMC3 modules, pass the `HEPMC3_PYTHON_VERSIONS` option to `cmake`, e.g.
```
cmake ... -DHEPMC3_PYTHON_VERSIONS=2,3.1,3.2
```
will build Python modules for versions Python 2.x (any Python of version 2), Python 3.1 and Python 3.2.
By default CMake will attempt to build the Python modules for both Python versions 2 and 3.


10b

In addition to the standard CPython  modules, it is possible to build HepMC3 modules for PyPy.
However, the PyPy support is experimental. To build the bindings against the `pypy-c` library use `pypy<version>`
for the `HEPMC3_PYTHON_VERSIONS` option, e.g.
```
-DHEPMC3_PYTHON_VERSIONS=pypy2
```

11

HepMC3 is shipped with interfaces to some MC event generators/codes located in the interfaces/ directory.
This is done to allow the usage of HepMC3 with codes that so far don't have HepMC3 interfaces.
In the future the codes from the interfaces  will be submitted to the
upstream of corresponding projects and removed from HepMC3.

To enable the installation of interfaces use `-DHEPMC3_INSTALL_INTERFACES:BOOL=ON`  option for the cmake.
In the test suite the already installed interfaces from MC event generators are used by default.
To use the interfaces and tests shipped with HepMC3, one can set the following options for cmake
```
-DHEPMC3_USE_INTERFACE_FROM_MCTESTER:BOOL=OFF
-DHEPMC3_USE_INTERFACE_FROM_PHOTOSPP:BOOL=OFF
-DHEPMC3_USE_INTERFACE_FROM_TAUOLAPP:BOOL=OFF
-DHEPMC3_USE_INTERFACE_FROM_PYTHIA8:BOOL=OFF
```

12.

 In order to build with PROTOBUFIO put the following flags
```
-DHEPMC3_ENABLE_PROTOBUF=ON 
```
Adjustment of LD_LIBRARY_PATH might be needed.


# Installation troubleshooting

The possible problems during the HepMC3 installation can be caused by

 - The used C++ compiler that does not support C++11.
   The only solution is to use compiler with C++11 support
 - The used CMake is too old.
 - Compilation of python bindings fails.
   pybind11 used by HepMC3 to create the python bindings supports only some compilers (gcc,clang, MSVC), therefore
   it is expected that the bindings will work  with these compilers.


#  Examples
HepMC3 is shipped with multiple examples. These can be compiled during installation as described above
or after the installation (for HepMC3>3.1.0).
To compile the examples after the installation  copy the installed directory with examples to desired directory
and run CMake, e.g.

```
mkdir -p myexamples
cd myexamples
cp -r /usr/share/doc/HepMC3-3.2.6/examples ./
cd examples
cmake -DUSE_INSTALLED_HEPMC3=ON CMakeLists.txt
make
```
Please note that in case the HepMC3 installation is not system-wide,
CMake will require an additional flag  `-DHepMC3_DIR=/where/you/have/it` to point to the HepMC3 installation directory.


# Source package structure

The package sources are organized as follows:

0. The topmost directory contains the
README and AUTHORS files, short design description, the main CMake file `CMakeLists.txt`, the package license,
the change log, and a template configuration file.
```
./DESIGN
./README-for-developers
./COPYING
./AUTHORS
./LICENCE
./HepMC3-config.in
./README.md
./ChangeLog
```

1. The subdirectory `src/` contains the core set of library sources while
the corresponding headers are located in include/HepMC3 directory.
The subdirectory `include/HepMC3/Data` includes headers with definitions of
POD structures used in the library.

2. The subdirectory `search/src` contains a set of source files of the search
sublibrary while the corresponding headers are located in
the `search/include/HepMC3` directory. The `search/CMakeLists.txt` file is used by CMake
to build  the search sublibrary.

3. The subdirectory `rootIO/src` contains a set of source files of the rootIO
sublibrary while the corresponding headers are located in
the `search/include/HepMC3` and `search/include/` directories. The
`rootIO/CMakeLists.txt` file is used by cmake to build  the search sublibrary.

4. The subdirectory `python` contains a set of source files required for compilation
of the Python bindings to the library. These include the `python/include/pybind11`
directory with the pybind11 headers.  The `python/src/` directory contains custom

```
python/src/root_binders.hpp
python/src/binders.hpp
python/src/stl_binders.hpp
```

and automatically generated binding sources

```
python/src/pyHepMC3*.cpp
python/src/pyHepMC3rootIO*.cpp
python/src/pyHepMC3protobufIO*.cpp
python/src/pyHepMC3search*.cpp
```

The automatically sources were generated using the binder configurations

```
python/src/pyHepMC3.binder
python/src/pyHepMC3rootIO.binder
python/src/pyHepMC3protobufIO.binder
python/src/pyHepMC3search.binder
```

and the header files generated from the templates

```
python/all_includes.hpp.in
python/root_includes.hpp.in
python/search_includes.hpp.in
```

The directory contain the files that are used to build a Python package and installation:
```
python/src/__init__.py
python/src/search/__init__.py
python/src/rootIO/__init__.py
python/src/protobufIO/__init__.py
pyHepMC3.egg-info.in
pyHepMC3.rootIO.egg-info.in
pyHepMC3.protobufIO.egg-info.in
pyHepMC3.search.egg-info.in
```

The `python/CMakeLists.txt` file is used by CMake
to build the binding.


5. The subdirectory `interfaces/` contains subdirectories with interfaces
(C++/Fortran source files, C++ headers) to
Monte Carlo event generators and event processing tools, e.g.
mc-tester,  Photospp,  pythia6,  pythia8 and  Tauolapp.


6. The subdirectories `test/` and `python/test` contain sets of files (source files, scripts, inputs)
that are used in the unit tests of the library and the Python bindings respectively.
The files in `test/` can be split in two groups:
- The tests of HepMC3 on itself with its inputs, e.g.
```
./test/testIO4.cc
./test/inputIO4.root
```
and those that are involved in testing of HepMC3 with Monte Carlo generators/tools

```
./test/McTesterValidationTool.cc
./test/PhotosValidationTool.h
./test/SimpleEventTool.cc
./test/McTesterValidationTool.h
./test/HepMC3TestUtils.h
./test/ValidationControl.cc
./test/OutputValidationTool.cc
./test/ValidationTool.h
./test/PythiaValidationTool.h
./test/TauolaValidationTool.cc
./test/SimpleEventTool.h
./test/PythiaValidationTool.cc
./test/OutputValidationTool.h
./test/IsGoodEvent.h
./test/ValidationControl.h
./test/Timer.h
./test/TauolaValidationTool.h
./test/PhotosValidationTool.cc
./test/testMCTester1.cc
./test/testPhotos1.cc
./test/testPythia3.cc
./test/testPythia2.cc
./test/testPythia1.cc
./test/testTauola1.cc
```


7. The subdirectory `cmake/Modules` contains files needed for the CMake configuration
The subdirectory `cmake/Templates` contains templates needed for generation
of the library CMake configuration files.


8. The subdirectory `examples` contains several examples of using the HepMC3 library
in standalone applications.
Each example is located in its own directory and can be built using the `CMakeLists.txt`
file in the same directory.

- The `examples/ConvertExample/` subdirectory contains source code of
utility that converts different types of event records, including cusm event records, into each other.
Namely, the files
```
./examples/ConvertExample/src/WriterDOT.cc
./examples/ConvertExample/include/WriterDOT.h
```
are the source files for an event format that can be visualized with graphviz.
The files
```
./examples/ConvertExample/src/WriterHEPEVTZEUS.cc
./examples/ConvertExample/include/WriterHEPEVTZEUS.h
```
contain an implementation of output format that can be used in the ZEUS experiment.

The files
```
./examples/ConvertExample/src/WriterRootTreeOPAL.cc
./examples/ConvertExample/include/WriterRootTreeOPAL.h
```
contain an implementation of output format that can be used together with data from the OPAL experiment.

The files
```
./examples/ConvertExample/src/AnalysisExample.cc
./examples/ConvertExample/include/AnalysisExample.h
```
illustrate an implementation of simple physics analysis using the HepMC3 library.

The files
```
./examples/ConvertExample/include/ReaderuprootTree.h
./examples/ConvertExample/src/ReaderuprootTree.cc
```
implement an a reader for ROOT files based on uproot.

- The `./examples/RootIOExample/` subdirectory contains source code of
an utility that illustrates manipulations with LHEF event record.

- The `./examples/RootIOExample/` subdirectory contains source code of
an utility that reads HepMC3 events in ROOT format.

- The `./examples/RootIOExample2/` subdirectory contains source code of
an utility that reads HepMC3 events in ROOT TTree format.

- The `./examples/RootIOExample2/` subdirectory contains source code of
an utility that reads HepMC3 events and saves them using a custom ROOT-based class.

- The `./examples/RootIOExample3/` subdirectory contains source code of
an utility that reads HepMC3 events and saves them nto ROOT TTree. This is a simplified version of ConvertExample.

- The `./examples/BasicExamples/` subdirectory contains source code of basic examples of HepMC3 usage,
e.g. building of event from scratch, reading and writing files, usage of fortran, etc.

- The `./examples/Pythia6Example/` subdirectory contains source code of
an utility that generates HepMC events with the Pythia6 Monte Carlo event
generator.

- The `./examples/Pythia8Example/` subdirectory contains source code of
an utility that generates HepMC events with the Pythia8 Monte Carlo event
generator.

- The `./examples/LHEFExample/` subdirectory contains source code of
an utility that illustrates manipulations with LHEF event record.

- The `./examples/ViewerExample/` subdirectory contains source code of ROOT
based GUI program that allows to visualize the HepMC3 events.

- The `./examples/SearchExample/` subdirectory contains source code example that
deals with search of relations between particles in the event.


9. The `doc/` directory contains files used for generation of library source code documentation with
the Doxygen system. The `doc/CMakeLists.txt` file is used by CMake
to build the documentation.
