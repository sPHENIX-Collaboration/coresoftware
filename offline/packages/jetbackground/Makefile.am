AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libjetbackground_io.la \
   libjetbackground.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libjetbackground_io_la_LIBADD = \
  -lphool

libjetbackground_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libjetbackground_la_LIBADD = \
  libjetbackground_io.la \
  -lcalo_io \
  -lConstituentSubtractor \
  -lg4jets_io \
  -lphg4hit \
  -lSubsysReco

pkginclude_HEADERS = \
  DetermineTowerBackground.h \
  SubtractTowers.h \
  SubtractTowersCS.h \
  RetowerCEMC.h \
  CopyAndSubtractJets.h \
  FastJetAlgoSub.h \
  TowerBackground.h \
  TowerBackgroundv1.h

ROOTDICTS = \
  TowerBackground_Dict.cc \
  TowerBackgroundv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TowerBackground_Dict_rdict.pcm \
  TowerBackgroundv1_Dict_rdict.pcm

libjetbackground_io_la_SOURCES = \
  $(ROOTDICTS) \
  TowerBackgroundv1.cc

libjetbackground_la_SOURCES = \
  CopyAndSubtractJets.cc \
  DetermineTowerBackground.cc \
  FastJetAlgoSub.cc \
  RetowerCEMC.cc \
  SubtractTowers.cc \
  SubtractTowersCS.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_jetbackground_io \
  testexternals_jetbackground


testexternals_jetbackground_io_SOURCES = testexternals.cc
testexternals_jetbackground_io_LDADD = libjetbackground_io.la

testexternals_jetbackground_SOURCES = testexternals.cc
testexternals_jetbackground_LDADD = libjetbackground.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
