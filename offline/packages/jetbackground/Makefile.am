AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

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
  -ljetbase \
  -lphg4hit \
  -lqautils \
  -lSubsysReco

pkginclude_HEADERS = \
  DetermineTowerBackground.h \
  DetermineTowerRho.h \
  SubtractTowers.h \
  SubtractTowersCS.h \
  RetowerCEMC.h \
  CopyAndSubtractJets.h \
  FastJetAlgoSub.h \
  TowerBackground.h \
  TowerBackgroundv1.h \
  TowerRho.h \
  TowerRhov1.h \
  BeamBackgroundFilterAndQA.h \
  BeamBackgroundFilterAndQADefs.h \
  BaseBeamBackgroundFilter.h \
  NullFilter.h \
  StreakSidebandFilter.h

ROOTDICTS = \
  TowerBackground_Dict.cc \
  TowerBackgroundv1_Dict.cc \
  TowerRho_Dict.cc \
  TowerRhov1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TowerBackground_Dict_rdict.pcm \
  TowerBackgroundv1_Dict_rdict.pcm \
  TowerRho_Dict_rdict.pcm \
  TowerRhov1_Dict_rdict.pcm

libjetbackground_io_la_SOURCES = \
  $(ROOTDICTS) \
  TowerBackgroundv1.cc \
  TowerRhov1.cc 

libjetbackground_la_SOURCES = \
  CopyAndSubtractJets.cc \
  DetermineTowerBackground.cc \
  DetermineTowerRho.cc \
  FastJetAlgoSub.cc \
  RetowerCEMC.cc \
  SubtractTowers.cc \
  SubtractTowersCS.cc \
  BeamBackgroundFilterAndQA.cc \
  NullFilter.cc \
  StreakSidebandFilter.cc

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
