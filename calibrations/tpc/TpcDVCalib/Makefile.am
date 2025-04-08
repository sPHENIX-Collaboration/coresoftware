##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  TrackToCalo.h

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libTpcDVCalib.la

libTpcDVCalib_la_SOURCES = \
  TrackToCalo.cc

libTpcDVCalib_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib

libTpcDVCalib_la_LIBADD = \
  -lActsCore \
  -lActsPluginTGeo \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework \
  -lSubsysReco \
  -ltrackbase_historic_io \
  -ltrack_io \
  -ltrack_reco \
  -lcalo_io \
  -lglobalvertex_io \
  -lg4detectors_io \
  -lphparameter_io

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES  = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libTpcDVCalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
