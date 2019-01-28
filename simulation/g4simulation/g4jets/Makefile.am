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
# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOTDICTS = \
  Jet_Dict.cc \
  JetV1_Dict.cc \
  JetMap_Dict.cc \
  JetMapV1_Dict.cc
# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  Jet_Dict_rdict.pcm \
  JetV1_Dict_rdict.pcm \
  JetMap_Dict_rdict.pcm \
  JetMapV1_Dict_rdict.pcm
else
  ROOT5DICTS = \
    ClusterJetInput_Dict.cc \
    FastJetAlgo_Dict.cc \
    JetReco_Dict.cc \
    JetHepMCLoader_Dict.cc \
    TowerJetInput_Dict.cc \
    TrackJetInput_Dict.cc \
    TruthJetInput_Dict.cc
endif

libg4jets_io_la_SOURCES = \
  $(ROOTDICTS) \
  Jet.cc \
  JetV1.cc \
  JetMap.cc \
  JetMapV1.cc

libg4jets_la_SOURCES = \
  $(ROOT5DICTS) \
  ClusterJetInput.cc \
  FastJetAlgo.cc \
  JetReco.cc \
  JetHepMCLoader.cc \
  TowerJetInput.cc \
  TrackJetInput.cc \
  TruthJetInput.cc

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
