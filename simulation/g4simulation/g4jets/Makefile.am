AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4jets_io.la \
   libg4jets.la

AM_CXXFLAGS = \
 -Wall -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libg4jets_io_la_LIBADD = \
  -lphool

libg4jets_la_LIBADD = \
  -lfun4all \
  -lphg4hit \
  -lg4hough_io \
  -lcalo_io -lcalo_util \
  -lg4vertex_io \
  -lCGAL \
  -lfastjet \
  -lphhepmc_io \
  -ltrackbase_historic_io \
  libg4jets_io.la

pkginclude_HEADERS = \
  Jet.h \
  JetV1.h \
  JetMap.h \
  JetMapV1.h \
  JetInput.h \
  JetAlgo.h \
  JetReco.h \
  TruthJetInput.h \
  TrackJetInput.h \
  TowerJetInput.h \
  ClusterJetInput.h \
  FastJetAlgo.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4jets_io_la_SOURCES = \
  Jet.C \
  Jet_Dict.C \
  JetV1.C \
  JetV1_Dict.C \
  JetMap.C \
  JetMap_Dict.C \
  JetMapV1.C \
  JetMapV1_Dict.C

libg4jets_la_SOURCES = \
  JetReco.C \
  JetReco_Dict.C \
  JetHepMCLoader.C \
  JetHepMCLoader_Dict.C \
  TruthJetInput.C \
  TruthJetInput_Dict.C \
  TrackJetInput.C \
  TrackJetInput_Dict.C \
  TowerJetInput.C \
  TowerJetInput_Dict.C \
  ClusterJetInput.C \
  ClusterJetInput_Dict.C \
  FastJetAlgo.C \
  FastJetAlgo_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4jets_io \
  testexternals_g4jets

testexternals_g4jets_io_SOURCES = testexternals.C
testexternals_g4jets_io_LDADD = libg4jets_io.la

testexternals_g4jets_SOURCES = testexternals.C
testexternals_g4jets_LDADD = libg4jets.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* testexternals.C *.pcm
