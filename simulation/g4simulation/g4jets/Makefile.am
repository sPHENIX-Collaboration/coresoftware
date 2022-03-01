AUTOMAKE_OPTIONS = foreign


AM_CXXFLAGS = `geant4-config --cflags`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4jets_io.la \
   libg4jets.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4jets_io_la_LIBADD = \
  -lphool

libg4jets_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libg4jets_la_LIBADD = \
  libg4jets_io.la \
  -lfun4all \
  -lphg4hit \
  -lcalo_io \
  -lg4vertex_io \
  -lRecursiveTools \
  -lphhepmc_io \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  Jet.h \
  Jetv1.h \
  JetMap.h \
  JetMapv1.h \
  JetInput.h \
  JetAlgo.h \
  JetReco.h \
  TruthJetInput.h \
  TrackJetInput.h \
  TowerJetInput.h \
  ClusterJetInput.h \
  FastJetAlgo.h \
  JetHepMCLoader.h

ROOTDICTS = \
  Jet_Dict.cc \
  Jetv1_Dict.cc \
  JetMap_Dict.cc \
  JetMapv1_Dict.cc

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
  Jet_Dict_rdict.pcm \
  Jetv1_Dict_rdict.pcm \
  JetMap_Dict_rdict.pcm \
  JetMapv1_Dict_rdict.pcm

libg4jets_io_la_SOURCES = \
  $(ROOTDICTS) \
  Jet.cc \
  Jetv1.cc \
  JetMap.cc \
  JetMapv1.cc

libg4jets_la_SOURCES = \
  ClusterJetInput.cc \
  FastJetAlgo.cc \
  JetReco.cc \
  JetHepMCLoader.cc \
  TowerJetInput.cc \
  TrackJetInput.cc \
  TruthJetInput.cc

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4jets_io \
  testexternals_g4jets

BUILT_SOURCES = testexternals.cc

testexternals_g4jets_io_SOURCES = testexternals.cc
testexternals_g4jets_io_LDADD = libg4jets_io.la

testexternals_g4jets_SOURCES = testexternals.cc
testexternals_g4jets_LDADD = libg4jets.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
