AUTOMAKE_OPTIONS = foreign

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(G4_MAIN)/include \
    -I$(G4_MAIN)/include/Geant4 \
    -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
    libTPCsimu.la

AM_CPPFLAGS = -DCGAL_DISABLE_ROUNDING_MATH_CHECK

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

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
  TPCCloud.h \
  TPCDetector.h \
  TPCSteppingAction.h \
  TPCEventAction.h \
  TPCDetectorSubsystem.h \
  TPCSimulation.h \
  TPCSimulationSubsystem.h 

libTPCsimu_la_SOURCES = \
  TPCCloud.cxx \
  TPCDetector.cxx \
  TPCSteppingAction.cxx \
  TPCEventAction.cxx \
  TPCDetectorSubsystem.cxx \
  TPCSimulation.cxx \
  TPCSimulationSubsystem.cxx \
  TPCsimu_Dict.cc

# Rule for generating table CINT dictionaries.
TPCsimu_Dict.cc: TPCDetectorSubsystem.h TPCSimulationSubsystem.h TPCsimuLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_TPCsimu
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
