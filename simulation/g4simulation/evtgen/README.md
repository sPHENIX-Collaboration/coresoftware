# EvtGen README

## Introduction
NEED TO WRITE SOMETHING

## Building the code
To simplify the procedure you can use the setupEvtGen.sh script, which
automates the process of downloading and building EvtGen and all dependencies
from source. There are variables in the script that can be edited to set the
installation location and the versions of EvtGen and each dependency (for more
details see below) to be used.

Alternatively you can perform a manual build, for example if you are building
against an existing installation of the dependencies, e.g. LCG releases.
This procedure is described in the rest of this section.

To build the EvtGen code, first make sure that there is a valid (C++) version
of HepMC avilable:

HepMC    http://hepmc.web.cern.ch/hepmc/

HepMC is used to store particle information.
We now recommend using HepMC3 but support for HepMC2 continues to be available.

Optionally, it is possible to use other external generators, such as Pythia8
(for Pythia decays in the DECAY.DEC file, for example), Photos (for
radiative corrections) and Tauola (for tau decays):

Pythia8  https://pythia.org/

Photos   http://photospp.web.cern.ch/photospp/

Tauola   http://tauolapp.web.cern.ch/tauolapp/

All of these packages have instructions for building them.

For HepMC3 support the following versions are required:
Pythia8: 8.201 or newer
Photos: 3.64 or newer
Tauola: 1.1.8 or newer

Once these packages are available, build the EvtGen release by creating a build directory
alongside the EvtGen source directory (assumed here to be called evtgen.git) and running:
```
cmake ../evtgen.git <options>
```
within the EvtGen build directory, using the following options:

  `-DCMAKE_INSTALL_PREFIX=<location>` : Location in which to install EvtGen (highly recommended)

  `-DEVTGEN_HEPMC3=ON`                : Enable HepMC3 support (default)
                                        To use HepMC2, set this to OFF instead.

  `-DHEPMC3_ROOT_DIR=<location>`      : Location of HepMC3 install directory
                                        While linking with either HepMC2 or HepMC3 is mandatory,
                                        depending on your environment the installation may be
                                        detected automatically.  Failing this, you can specify
                                        the location via this option.

  `-DHEPMC2_ROOT_DIR=<location>`      : Location of HepMC2 install directory
                                        While linking with either HepMC2 or HepMC3 is mandatory,
                                        depending on your environment the installation may be
                                        detected automatically.  Failing this, you can specify
                                        the location via this option.

  `-DEVTGEN_PYTHIA=ON`                : Enable linking with Pythia 8 (OFF by default)

  `-DPYTHIA8_ROOT_DIR=<location>`     : Location of Pythia8 install directory
                                        As with HepMC this may be automatically detected
                                        depending on your build environment, otherwise the
                                        location can be specified via this option.

  `-DEVTGEN_PHOTOS=ON`                : Enable linking with Photos++ (OFF by default)

  `-DPhotos++_ROOT_DIR=<location>`
  or
  `-DPHOTOSPP_ROOT_DIR=<location>`    : Location of Photos++ install directory
                                        As with HepMC this may be automatically detected
                                        depending on your build environment, otherwise the
                                        location can be specified via this option.

  `-DEVTGEN_TAUOLA=ON`                : Enable linking with Tauola++ (OFF by default)

  `-DTauola++_ROOT_DIR=<location>`
  or
  `-DTAUOLAPP_ROOT_DIR=<location>`    : Location of Tauola++ install directory
                                        As with HepMC this may be automatically detected
                                        depending on your build environment, otherwise the
                                        location can be specified via this option.

  `-DEVTGEN_BUILD_DOC=ON`             : Enable building documentation in 'doc' directory (OFF by default)

  `-DEVTGEN_BUILD_TESTS=ON`           : Enable building executables in 'test' directory (OFF by default)

  `-DEVTGEN_BUILD_VALIDATIONS=ON`     : Enable building executables in 'validation' directory (OFF by default)

Then compile and (optionally, although highly recommended) install the EvtGen code using
```
make
make install
```
This should create the libraries libEvtGen.so and libEvtGenExternal.so, as well as the archives
libEvtGen.a and libEvtGenExternal.a (the "EvtGenExternal" library/archive will not be created
if linking is not enabled for any of the external generators).

A series of validation and test executables are also created (discussed below).
The CMakeLists.txt file in these two directories can be used as a template for building any
user executables.

The build/install will also create CMake config files that allow the library targets to be
easily imported into other projects.

To use Pythia 8, the environment variable PYTHIA8DATA needs to be set to the location
of the corresponding xml documentation directory, which also contains the default values
for particle decays and models: `<location of Pythia 8 base directory>/xmldoc`


## Code validation
The validation sub-directory contains code validation test cases for Pythia, Tauola
and B mixing models.

Note that these executables also depend on ROOT, which is only used to create ntuples and
plots for the validation examples.

The script genAllDecayExamples.sh runs other scripts that generate a range of decay
modes using the genExampleRootFiles.cc program. The script compareAllDecays.sh runs
a range of scripts that creates comparison plots using the compareRootFiles.cc program.
For now, the comparisons use the same plots, but each of the "compare.sh" files can
be edited to compare any two versions of ROOT data created by the
genExampleRootFiles.cc program.

The testCPVDecays.cc program runs a test for the B mixing decay model.


## Examples
Some examples are provided in the test sub-directory.

Note that these executables also depend on ROOT, which is only used to create ntuples and
plots for the validation examples.

Running the script `./do_tests` will run a series of EvtGen examples.
Example decay files are in the test/exampleFiles sub-directory.


## Release notes
Please see [History.md](History.md) for a detailed list of changes to this package.

The major points comparing this version with the 2009 release are the following:

1. This version requires HepMC (version 2.04 and above) for storing event structures
   of particle decays, and can also use Pythia 8 (version 8.180 and above is required)
   and the C++ interfaced packages Photos (version 3.5.2 and above) and Tauola
   (version 1.0.7 and above). These external generators are included via engine
   classes in the new sub-directory EvtGenExternal.

   Two libraries can be created for EvtGen:
   i) libEvtGen.so contains the EvtGenBase and EvtGenModel core code/decay models,
   ii) libEvtGenExternal.so contains _only_ the code within EvtGenExternal.
   This means that the external generator interface can be ignored by not
   loading/creating the 2nd library libEvtGenExternal.so.

   In the installation instructions below, it is possible to select which external
   generators you want to use. Note, however, that the generic "DECAY.DEC" file
   contains Pythia decays, and if Pythia is not included in the build, these
   decays are not generated (in fact, they would need a new decay model specified).
   Likewise, if Photos is not included, there will be no radiative corrections
   done for particle decays.

   To use the external generators, use the following code:
   ```c++
   #include "EvtGenExternal/EvtExternalGenList.hh"
   #include "EvtGenBase/EvtAbsRadCorr.hh"
   #include "EvtGenBase/EvtDecayBase.hh"

   // Set up the default external generator list: Photos, Pythia and/or Tauola
   EvtExternalGenList genList;
   EvtAbsRadCorr* radCorrEngine = genList.getPhotosModel();
   std::list<EvtDecayBase*> extraModels = genList.getListOfModels();

   // Create the EvtGen generator object
   EvtGen myGenerator("decayFile.dec", "evt.pdl", randomEnginePointer,
                      radCorrEngine, &extraModels);

   //If you don't want to use external generators, use the following:
   //EvtGen myGenerator("decayFile.dec", "evt.pdl", randomEnginePointer);
   ```

   The files [Pythia8_README.md](Pythia8_README.md) and [Tauola_README.md](Tauola_README.md) have more details about
   using the new Pythia 8 and Tauola generators (called via the PYTHIA and TAUOLA
   "decay.dec" model names). The new Photos generator is still called via the PHOTOS
   "decay.dec" model name.

   It is now possible to use alias particle decays for the Pythia 8 model.
   Two Pythia 8 instances are used in EvtPythiaEngine for normal and aliased decays.
   Since the underlying code for Photos and Tauola is still Fortran, it is only
   possible to have one (unique) instance of each of these external generators.
   This can only be fixed if these packages are converted to pure C++ code.


2. This version of EvtGen is effectively a merger of the latest LHCb and BaBar
   EvtGenBase and EvtGenModels code. Various decay models have been added, and
   there have been a range of bug fixes.


3. There is also a new Dalitz decay class model "GENERIC_DALITZ"
   (EvtGenModels/EvtGenericDalitz.cpp), that should be used instead of EvtDDalitz.
   The generic Dalitz model uses xml files to configure the resonance amplitude
   parameters (instead of being hardcoded in EvtDDalitz):

   ```
   Decay D+
   1.0  K-  pi+  pi+  GENERIC_DALITZ MyDalitzParameters.xml;
   Enddecay
   ```

   Examples of xml Dalitz parameter files are given in the sub-directory
   validation/DalitzFiles, e.g. see DalitzDecays.xml.


4. It is possible to use decay files in xml format. Use the python script
   convertDecayFile.py for converting decay files to the new format.
   The src/EvtGen.cpp constructor has an additional boolean argument useXml that
   needs to be set to true if xml decay files are to be used (default is useXml=false).
   For example, DECAY_2010.XML is the xml version of DECAY_2010.DEC.


5. Bug fixes for Bs mixing decay/CP violation amplitudes.
   Added the capability to use either coherent or incoherent mixing in EvtCPUtil.
   One or the other can be chosen as the mixing method for the B system by
   choosing 0 (coherent) or 1 (incoherent) for the last integer argument
   in the EvtGen() constructor.


