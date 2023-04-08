##############################################
# please add new classes/includes/everything else in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libTPCHitTrackDisplay.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(G4_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libTPCHitTrackDisplay_la_LIBADD = \
  -lphool \
  -lg4detectors \
  -lcentrality_io \
  -lphg4hit \
  -lConstituentSubtractor \
  -lg4jets_io \
  -ltrackbase_historic_io \
  -lSubsysReco \
  -lcalo_io

pkginclude_HEADERS = \
  TPCHitTrackDisplay.h

libTPCHitTrackDisplay_la_SOURCES = \
  TPCHitTrackDisplay.cc

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libTPCHitTrackDisplay.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
