AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libjetbase_io.la \
   libjetbase.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib 

libjetbase_io_la_LIBADD = \
  -lphool

libjetbase_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs` \
  -lConstituentSubtractor

libjetbase_la_LIBADD = \
  libjetbase_io.la \
  -lcalo_io \
  -lcdbobjects \
  -lffamodules \
  -lfun4all \
  -lglobalvertex_io \
  -lgsl \
  -lgslcblas \
  -lRecursiveTools \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  ClusterJetInput.h \
  FastJetAlgo.h \
  FastJetOptions.h \
  Jet.h \
  JetCalib.h \
  Jetv1.h \
  Jetv2.h \
  JetContainer.h \
  JetContainerv1.h \
  JetMap.h \
  JetMapv1.h \
  JetInput.h \
  JetProbeMaker.h \
  JetProbeInput.h \
  JetAlgo.h \
  JetReco.h \
  TowerJetInput.h \
  TrackJetInput.h

ROOTDICTS = \
  Jet_Dict.cc \
  Jetv1_Dict.cc \
  Jetv2_Dict.cc \
  JetContainer_Dict.cc \
  JetContainerv1_Dict.cc \
  JetMap_Dict.cc \
  JetMapv1_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libjetbase_io_la_SOURCES = \
  $(ROOTDICTS) \
  Jet.cc \
  Jetv1.cc \
  Jetv2.cc \
  JetContainer.cc \
  JetContainerv1.cc \
  JetMap.cc \
  JetMapv1.cc

libjetbase_la_SOURCES = \
  ClusterJetInput.cc \
  JetAlgo.cc \
  FastJetAlgo.cc \
  FastJetOptions.cc \
  JetCalib.cc \
  JetProbeMaker.cc \
  JetProbeInput.cc \
  JetReco.cc \
  TowerJetInput.cc \
  TrackJetInput.cc

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_jetbase_io \
  testexternals_jetbase

BUILT_SOURCES = testexternals.cc

testexternals_jetbase_io_SOURCES = testexternals.cc
testexternals_jetbase_io_LDADD = libjetbase_io.la

testexternals_jetbase_SOURCES = testexternals.cc
testexternals_jetbase_LDADD = libjetbase.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
