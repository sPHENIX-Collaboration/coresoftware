##############################################
# please add new classes/includes/everything else in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libg4tpc.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4tpc_la_LIBADD = \
  -lphool \
  -lg4detectors \
  -lphg4hit \
  -lphparameter \
  -ltpc_io

pkginclude_HEADERS = \
  PHG4TpcCentralMembrane.h \
  TpcClusterBuilder.h \
  PHG4TpcDigitizer.h \
  PHG4TpcDirectLaser.h \
  PHG4TpcDistortion.h \
  PHG4TpcElectronDrift.h \
  PHG4TpcEndCapSubsystem.h \
  PHG4TpcPadBaselineShift.h \
  PHG4TpcPadPlane.h \
  PHG4TpcPadPlaneReadout.h \
  PHG4TpcSubsystem.h 

libg4tpc_la_SOURCES = \
  PHG4TpcCentralMembrane.cc \
  TpcClusterBuilder.cc \
  PHG4TpcDetector.cc \
  PHG4TpcDigitizer.cc \
  PHG4TpcDirectLaser.cc \
  PHG4TpcDisplayAction.cc \
  PHG4TpcDistortion.cc \
  PHG4TpcElectronDrift.cc \
  PHG4TpcEndCapDetector.cc \
  PHG4TpcEndCapDisplayAction.cc \
  PHG4TpcEndCapSteppingAction.cc \
  PHG4TpcEndCapSubsystem.cc \
  PHG4TpcPadBaselineShift.cc \
  PHG4TpcPadPlane.cc \
  PHG4TpcPadPlaneReadout.cc \
  PHG4TpcSteppingAction.cc \
  PHG4TpcSubsystem.cc 

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libg4tpc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
