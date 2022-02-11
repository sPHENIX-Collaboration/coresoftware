AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
  libg4trackfastsim.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4trackfastsim_la_LIBADD = \
  -lcalo_io \
  -lgsl \
  -lgslcblas \
  -lPHGenFit \
  -lphgeom \
  -lphg4hit \
  -ltrackbase_historic_io \
  -lphparameter

pkginclude_HEADERS = \
  PHG4TrackFastSim.h \
  PHG4TrackFastSimEval.h

libg4trackfastsim_la_SOURCES = \
  PHG4TrackFastSim.cc \
  PHG4TrackFastSimEval.cc


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.cc
testexternals_tracker_LDADD = libg4trackfastsim.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f $(BUILT_SOURCES)
