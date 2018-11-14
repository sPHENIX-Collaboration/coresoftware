AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libg4histos.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I$(G4_MAIN)/include

noinst_HEADERS = \
  G4AnalysisLinkDef.h \
  G4CellNtuple.h \
  G4EdepNtuple.h \
  G4EvtTree.h \
  G4HitTTree.h \
  G4RawTowerTTree.h \
  G4RootHitContainer.h \
  G4RootRawTower.h \
  G4RootRawTowerContainer.h \
  G4SnglNtuple.h \
  G4SnglTree.h

pkginclude_HEADERS = \
  G4HitNtuple.h

libg4histos_la_SOURCES = \
  G4CellNtuple.cc \
  G4EdepNtuple.cc \
  G4HitNtuple.cc \
  G4HitTTree.cc \
  G4RawTowerTTree.cc \
  G4RootHitContainer.cc \
  G4RootRawTower.cc \
  G4RootRawTowerContainer.cc \
  G4RootScintillatorSlat.cc \
  G4RootScintillatorSlatContainer.cc \
  G4ScintillatorSlatTTree.cc \
  G4RootScintillatorTower.cc \
  G4RootScintillatorTowerContainer.cc \
  G4ScintillatorTowerTTree.cc \
  G4SnglNtuple.cc \
  G4SnglTree.cc \
  G4TowerNtuple.cc

nodist_libg4histos_la_SOURCES = \
  G4Analysis_Dict.C

libg4histos_la_LDFLAGS = \
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
testexternals_LDADD = libg4histos.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# rules to make root dictionaries
G4Analysis_Dict.C: \
  G4CellNtuple.h \
  G4EdepNtuple.h \
  G4EvtTree.h \
  G4HitNtuple.h \
  G4HitTTree.h \
  G4RawTowerTTree.h \
  G4RootHitContainer.h \
  G4RootRawTower.h \
  G4RootRawTowerContainer.h \
  G4RootScintillatorSlat.h \
  G4RootScintillatorSlatContainer.h \
  G4ScintillatorSlatTTree.h \
  G4RootScintillatorTower.h \
  G4RootScintillatorTowerContainer.h \
  G4ScintillatorTowerTTree.h \
  G4SnglNtuple.h \
  G4SnglTree.h \
  G4TowerNtuple.h \
  G4AnalysisLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
