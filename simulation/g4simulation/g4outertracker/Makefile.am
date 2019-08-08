AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libg4outertracker.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4outertracker_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4detectors \
  -ltrack_io \
  -loutertracker

pkginclude_HEADERS = \
  PHG4OuterTrackerDefs.h \
  PHG4OuterTrackerHitReco.h \
  PHG4OuterTrackerSubsystem.h \
  PHG4OuterTrackerDigitizer.h 

if MAKEROOT6
else
  ROOT5_DICTS = \
    PHG4OuterTrackerHitReco_Dict.cc \
    PHG4OuterTrackerDefs_Dict.cc \
    PHG4OuterTrackerSubsystem_Dict.cc \
    PHG4OuterTrackerDigitizer_Dict.cc
endif

libg4outertracker_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4OuterTrackerHitReco.cc \
  PHG4OuterTrackerDetector.cc \
  PHG4OuterTrackerDisplayAction.cc \
  PHG4OuterTrackerSteppingAction.cc \
  PHG4OuterTrackerSubsystem.cc \
  PHG4OuterTrackerDigitizer.cc

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

testexternals_SOURCES = \
  testexternals.cc
testexternals_LDADD = \
  libg4outertracker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
