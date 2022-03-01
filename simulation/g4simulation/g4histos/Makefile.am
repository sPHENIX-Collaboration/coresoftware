AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

lib_LTLIBRARIES = \
    libg4histos.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4histos_la_LDFLAGS = \
   -L$(libdir) \
   -L$(OFFLINE_MAIN)/lib \
   -lcalo_io \
   -lfun4all \
   -lg4detectors_io \
   -lphg4hit


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

noinst_HEADERS = \
  G4RootHitContainer.h \
  G4RootRawTower.h \
  G4RootRawTowerContainer.h

pkginclude_HEADERS = \
  G4CellNtuple.h \
  G4EdepNtuple.h \
  G4EvtTree.h \
  G4HitNtuple.h \
  G4HitTTree.h \
  G4RawTowerTTree.h \
  G4SnglNtuple.h \
  G4SnglTree.h \
  G4VtxNtuple.h

ROOTDICTS = \
  G4RootHitContainer_Dict.cc \
  G4RootRawTowerContainer_Dict.cc \
  G4RootRawTower_Dict.cc \
  G4RootScintillatorSlatContainer_Dict.cc \
  G4RootScintillatorSlat_Dict.cc \
  G4RootScintillatorTowerContainer_Dict.cc \
  G4RootScintillatorTower_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  G4RootHitContainer_Dict_rdict.pcm \
  G4RootRawTowerContainer_Dict_rdict.pcm \
  G4RootRawTower_Dict_rdict.pcm \
  G4RootScintillatorSlatContainer_Dict_rdict.pcm \
  G4RootScintillatorSlat_Dict_rdict.pcm \
  G4RootScintillatorTowerContainer_Dict_rdict.pcm \
  G4RootScintillatorTower_Dict_rdict.pcm

libg4histos_la_SOURCES = \
  $(ROOTDICTS) \
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
  G4TowerNtuple.cc \
  G4VtxNtuple.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libg4histos.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
