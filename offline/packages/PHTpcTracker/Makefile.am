##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libPHTpcTracker.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS =


libPHTpcTracker_la_SOURCES = \
  Fitter.cc \
  PHTpcEventExporter.cc \
  PHTpcLookup.cc \
  PHTpcSeedFinder.cc \
  PHTpcTrackFollower.cc \
  PHTpcTracker.cc \
  PHTpcTrackerUtil.cc \
  PHTpcVertexFinder.cc \
  SpacepointMeasurement2.cc \
  Track.cc


libPHTpcTracker_la_LIBADD = \
  -llog4cpp



# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals


testexternals_SOURCES = testexternals.cc
testexternals_track_reco_LDADD = libPHTpcTracker.la

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
