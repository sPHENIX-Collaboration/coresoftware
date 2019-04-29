AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4mvtx.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4mvtx_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors \
  -ltrack_io \
  -lmvtx_io

pkginclude_HEADERS = \
  PHG4MVTXHitReco.h \
  PHG4MVTXSubsystem.h \
  PHG4MvtxDigitizer.h 

if MAKEROOT6
else
  ROOT5_DICTS = \
    PHG4MVTXHitReco_Dict.cc \
    PHG4MVTXDefs_Dict.cc \
    PHG4MVTXSubsystem_Dict.cc \
    PHG4MvtxDigitizer_Dict.cc
endif

libg4mvtx_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4MVTXHitReco.cc \
  PHG4MVTXDetector.cc \
  PHG4MVTXSteppingAction.cc \
  PHG4MVTXSubsystem.cc \
  PHG4MvtxDigitizer.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc
testexternals_LDADD = \
  libg4mvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
