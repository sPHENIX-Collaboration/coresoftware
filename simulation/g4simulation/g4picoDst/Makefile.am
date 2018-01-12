AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libg4picoDst.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

noinst_HEADERS = \
  G4AnalysisLinkDef.h \
  G4SnglHit.h \
  G4SnglTower.h \
  G4SnglCluster.h \
  G4SnglHTCContainer.h \
  mFillG4SnglHTCContainer.h

libg4picoDst_la_SOURCES = \
  G4SnglHit.cc \
  G4SnglTower.cc \
  G4SnglCluster.cc \
  G4SnglHTCContainer.cc \
  mFillG4SnglHTCContainer.cc

nodist_libg4picoDst_la_SOURCES = \
  G4Analysis_Dict.C

libg4picoDst_la_LDFLAGS = \
   -L$(libdir) \
   -L$(OFFLINE_MAIN)/lib \
   -lcalo_util \
   -lfun4all \
   -lphg4hit \
   -lg4detectors \
   -lg4hough \
   -lg4testbench


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libg4picoDst.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

CLEANFILES = G4Analysis_Dict.h G4Analysis_Dict.C

# rules to make root dictionaries
G4Analysis_Dict.C: \
  G4SnglHit.h \
  G4SnglTower.h \
  G4SnglCluster.h \
  G4SnglHTCContainer.h \
  mFillG4SnglHTCContainer.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

