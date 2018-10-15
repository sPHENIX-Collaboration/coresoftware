AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libg4tpc_io.la \
    libg4tpc.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(G4_MAIN)/include \
    -I$(ROOTSYS)/include

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4tpc_la_LIBADD = \
  libg4tpc_io.la \
  -lphool \
  -lg4testbench \
  -lphg4hit \
  -lg4detectors \
  -lphparameter

pkginclude_HEADERS = \
  PHG4CellTPCv1.h

libg4tpc_la_SOURCES = \
  PHG4TPCDetector.cc \
  PHG4TPCElectronDrift.cc \
  PHG4TPCElectronDrift_Dict.cc \
  PHG4TPCPadPlane.cc \
  PHG4TPCPadPlane_Dict.cc \
  PHG4TPCPadPlaneSimple.cc \
  PHG4TPCPadPlaneSimple_Dict.cc \
  PHG4TPCPadPlaneReadout.cc \
  PHG4TPCPadPlaneReadout_Dict.cc \
  PHG4TPCSteppingAction.cc \
  PHG4TPCSubsystem.cc \
  PHG4TPCSubsystem_Dict.cc

libg4tpc_io_la_SOURCES = \
  PHG4CellTPCv1.cc \
  PHG4CellTPCv1_Dict.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


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

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
