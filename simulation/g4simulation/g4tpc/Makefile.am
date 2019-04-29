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
  PHG4TPCPadPlane.h \
  PHG4TPCPadPlaneReadout.h \
  PHG4TpcDigitizer.h \
  PHG4TPCSubsystem.h \
  PHG4TPCDistortion.h \
  PHG4TPCSpaceChargeDistortion.h 

if MAKEROOT6
else
  ROOT5_DICTS = \
    PHG4TpcElectronDrift_Dict.cc \
    PHG4TPCPadPlane_Dict.cc \
    PHG4TPCPadPlaneReadout_Dict.cc \
    PHG4TpcDigitizer_Dict.cc \
    PHG4TPCSubsystem_Dict.cc \
    PHG4TPCDistortion_Dict.cc\
    PHG4TPCSpaceChargeDistortion_Dict.cc 
endif

libg4tpc_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4TPCDetector.cc \
  PHG4TpcElectronDrift.cc \
  PHG4TPCPadPlane.cc \
  PHG4TPCPadPlaneReadout.cc \
  PHG4TPCSteppingAction.cc \
  PHG4TpcDigitizer.cc \
  PHG4TPCSubsystem.cc \
  PHG4TPCDistortion.cc \
  PHG4TPCSpaceChargeDistortion.cc

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

testexternals_io_SOURCES = testexternals.cc

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
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
