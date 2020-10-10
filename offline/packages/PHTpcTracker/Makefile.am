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

pkginclude_HEADERS = \
  PHTpcTracker.h


##############################################
# please add new classes in alphabetical order

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
  -lgenfit2 \
  -lgenfit2exp \
  -llog4cpp \
  -lphool \
  -lphfield \
  -lphgeom \
  -lSubsysReco \
  -ltrack_reco \
  -ltrack_io

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals


testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libPHTpcTracker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
