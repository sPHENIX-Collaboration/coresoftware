##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  liboutertracker_io.la \
  liboutertracker.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  OuterTrackerDefs.h \
  OuterTrackerHit.h \
  OuterTrackerClusterizer.h \
  CylinderGeomOuterTracker.h

ROOTDICTS = \
  OuterTrackerHit_Dict.cc \
  CylinderGeomOuterTracker_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  OuterTrackerHit_Dict_rdict.pcm \
  CylinderGeomOuterTracker_Dict_rdict.pcm
else
  ROOT5_DICTS = \
  OuterTrackerClusterizer_Dict.cc

endif

# sources for outertracker library
liboutertracker_la_SOURCES = \
  $(ROOT5_DICTS) \
  OuterTrackerClusterizer.cc

liboutertracker_la_LIBADD = \
  liboutertracker_io.la \
  -lCLHEP \
  -lphg4hit \
  -lSubsysReco

# sources for io library
liboutertracker_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeomOuterTracker.cc \
  OuterTrackerDefs.cc \
  OuterTrackerHit.cc

liboutertracker_io_la_LIBADD = \
  -lCLHEP \
  -lg4detectors_io \
  -lphool \
  -ltrack_io \
  -ltrackbase_historic_io


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_outertracker_io \
  testexternals_outertracker

testexternals_outertracker_io_SOURCES = testexternals.cc
testexternals_outertracker_io_LDADD = liboutertracker_io.la

testexternals_outertracker_SOURCES = testexternals.cc
testexternals_outertracker_LDADD = liboutertracker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
