##############################################
# please add new classes/includes/everything else in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
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
  -lphool \
  -lg4detectors \
  -lphg4hit \
  -lphparameter \
  -ltpc_io

pkginclude_HEADERS = \
  PHG4TpcElectronDrift.h \
  PHG4TpcPadPlane.h \
  PHG4TpcPadPlaneReadout.h \
  PHG4TpcDigitizer.h \
  PHG4TpcSubsystem.h \
  PHG4TpcDistortion.h \
  PHG4TpcSpaceChargeDistortion.h 

libg4tpc_la_SOURCES = \
  PHG4TpcDetector.cc \
  PHG4TpcDigitizer.cc \
  PHG4TpcDisplayAction.cc \
  PHG4TpcDistortion.cc \
  PHG4TpcElectronDrift.cc \
  PHG4TpcPadPlane.cc \
  PHG4TpcPadPlaneReadout.cc \
  PHG4TpcSpaceChargeDistortion.cc \
  PHG4TpcSteppingAction.cc \
  PHG4TpcSubsystem.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

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
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
