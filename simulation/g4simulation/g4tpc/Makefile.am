AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libg4tpc_io.la \
    libg4tpc.la

AM_CPPFLAGS = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I$(G4_MAIN)/include \
    -I$(ROOTSYS)/include \
    -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
    -L$(libdir) \
    -L$(OFFLINE_MAIN)/lib

libg4tpc_la_LIBADD = \
  libg4tpc_io.la \
  -lphool \
  -lg4testbench \
  -lphg4hit \
  -lg4detectors \
  -lphparameter \
  -ltrackbase_historic_io \
  -ltrack_io \
  -ltpc

pkginclude_HEADERS = \
  PHG4CellTPCv1.h \
  PHG4TPCElectronDrift.h \
  PHG4TPCPadPlane.h \
  PHG4TPCPadPlaneReadout.h \
  PHG4TPCDigitizer.h \
  PHG4TpcDigitizer.h \
  PHG4TPCSubsystem.h \
  PHG4TPCDistortion.h \
  PHG4TPCSpaceChargeDistortion.h 

ROOTDICTS = \
  PHG4CellTPCv1_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    PHG4CellTPCv1_Dict_rdict.pcm
else
  ROOT5_DICTS = \
    PHG4TPCElectronDrift_Dict.cc \
    PHG4TPCPadPlane_Dict.cc \
    PHG4TPCPadPlaneSimple_Dict.cc \
    PHG4TPCPadPlaneReadout_Dict.cc \
    PHG4TPCDigitizer_Dict.cc \
    PHG4TpcDigitizer_Dict.cc \
    PHG4TPCSubsystem_Dict.cc \
    PHG4TPCDistortion_Dict.cc\
    PHG4TPCSpaceChargeDistortion_Dict.cc 
endif

libg4tpc_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4TPCDetector.cc \
  PHG4TPCElectronDrift.cc \
  PHG4TPCPadPlane.cc \
  PHG4TPCPadPlaneSimple.cc \
  PHG4TPCPadPlaneReadout.cc \
  PHG4TPCSteppingAction.cc \
  PHG4TPCDigitizer.cc \
  PHG4TpcDigitizer.cc \
  PHG4TPCSubsystem.cc \
  PHG4TPCDistortion.cc \
  PHG4TPCSpaceChargeDistortion.cc

libg4tpc_io_la_SOURCES = \
  $(ROOTDICTS) \
  PHG4CellTPCv1.cc

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

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
