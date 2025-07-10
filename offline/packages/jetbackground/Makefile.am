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
  -leventplaneinfo_io \
  -lglobalvertex \
  -ljetbase \
  -lpdbcalBase \
  -lphg4hit \
  -lphparameter \
  -lqautils \
  -lSubsysReco

pkginclude_HEADERS = \
  BeamBackgroundFilterAndQA.h \
  BeamBackgroundFilterAndQADefs.h \
  BaseBeamBackgroundFilter.h \
  CopyAndSubtractJets.h \
  DetermineTowerBackground.h \
  DetermineEventRho.h \
  DetermineTowerRho.h \
  EventRho.h \
  EventRhov1.h \
  FastJetAlgoSub.h \
  JetBackgroundCut.h \
  NullFilter.h \
  RetowerCEMC.h \
  StreakSidebandFilter.h \
  SubtractTowers.h \
  SubtractTowersCS.h \
  TowerBackground.h \
  TowerBackgroundv1.h \
  TowerRho.h \
  TowerRhov1.h

ROOTDICTS = \
  EventRho_Dict.cc \
  EventRhov1_Dict.cc \
  TowerBackground_Dict.cc \
  TowerBackgroundv1_Dict.cc \
  TowerRho_Dict.cc \
  TowerRhov1_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libjetbackground_io_la_SOURCES = \
  $(ROOTDICTS) \
  EventRhov1.cc \
  TowerBackgroundv1.cc \
  TowerRhov1.cc

libjetbackground_la_SOURCES = \
  BeamBackgroundFilterAndQA.cc \
  CopyAndSubtractJets.cc \
  DetermineTowerBackground.cc \
  DetermineEventRho.cc \
  DetermineTowerRho.cc \
  FastJetAlgoSub.cc \
  JetBackgroundCut.cc \
  NullFilter.cc \
  RetowerCEMC.cc \
  StreakSidebandFilter.cc \
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
