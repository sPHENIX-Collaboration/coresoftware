AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(G4_MAIN)/include \
    -I$(G4_MAIN)/include/Geant4 \
    -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
    libTPCbase.la \
    libTPCsimu.la

AM_CPPFLAGS = -DCGAL_DISABLE_ROUNDING_MATH_CHECK

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libTPCbase_la_LIBADD = \
  -L$(ROOTSYS)/lib \
  -lphool  \
  -lboost_filesystem \
  -lboost_system \
  -lXMLIO \
  -lCGAL \
  -lphg4hit

libTPCsimu_la_LIBADD = \
  -L$(ROOTSYS)/lib \
  -lphool  \
  -lboost_filesystem \
  -lboost_system \
  -lpdbcalBase \
  -lXMLIO \
  -lCGAL \
  -lSubsysReco \
  -lg4testbench \
  -lphg4hit \
  -lg4detectors \
  -lTPCbase

pkginclude_HEADERS = \
  TPCDataTypes.h \
  TPCConstants.h \
  TPCPadMap.h \
  TPCHit.h \
  TPCDigit.h \
  TPCHitsContainer.h \
  TPCDigitsContainer.h \
  TPCElectronPDF.h \
  TPCDetector.h \
  TPCSteppingAction.h \
  TPCEventAction.h \
  TPCDetectorSubsystem.h \
  TPCSimulation.h \
  TPCSimulationSubsystem.h 

libTPCbase_la_SOURCES = \
  TPCPadMap.cxx \
  TPCHit.cxx \
  TPCDigit.cxx \
  TPCHitsContainer.cxx \
  TPCDigitsContainer.cxx \
  TPCbase_Dict.cc

libTPCsimu_la_SOURCES = \
  TPCElectronPDF.cxx \
  TPCDetector.cxx \
  TPCSteppingAction.cxx \
  TPCEventAction.cxx \
  TPCDetectorSubsystem.cxx \
  TPCSimulation.cxx \
  TPCSimulationSubsystem.cxx \
  TPCsimu_Dict.cc

# Rule for generating table CINT dictionaries.
TPCbase_Dict.cc: TPCPadMap.h TPCHit.h TPCDigit.h TPCHitsContainer.h TPCDigitsContainer.h TPCbaseLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
TPCsimu_Dict.cc: TPCDetectorSubsystem.h TPCSimulationSubsystem.h TPCsimuLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_TPCbase \
  testexternals_TPCsimu

testexternals_TPCbase_SOURCES = testexternals.cc
testexternals_TPCbase_LDADD = libTPCbase.la

testexternals_TPCsimu_SOURCES = testexternals.cc
testexternals_TPCsimu_LDADD = libTPCsimu.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
