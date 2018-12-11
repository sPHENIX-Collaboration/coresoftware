AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libjetbackground_io.la \
   libjetbackground.la

AM_CXXFLAGS = \
 -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs`

libjetbackground_io_la_LIBADD = \
  -lphool

libjetbackground_la_LIBADD = \
  -lfun4all \
  -lcalo_io -lcalo_util \
  -lg4jets \
  -lg4jets_io \
  -lCGAL \
  -lfastjet \
  -lfastjettools \
  -lConstituentSubtractor \
  libjetbackground_io.la

pkginclude_HEADERS = \
  DetermineTowerBackground.h \
  SubtractTowers.h \
  SubtractTowersCS.h \
  RetowerCEMC.h \
  CopyAndSubtractJets.h \
  FastJetAlgoSub.h \
  TowerBackground.h \
  TowerBackground_v1.h

libjetbackground_io_la_SOURCES = \
  TowerBackground_Dict.C \
  TowerBackground_v1.C \
  TowerBackground_v1_Dict.C

libjetbackground_la_SOURCES = \
  DetermineTowerBackground.C \
  DetermineTowerBackground_Dict.C \
  FastJetAlgoSub.C \
  FastJetAlgoSub_Dict.C \
  SubtractTowers.C \
  SubtractTowers_Dict.C \
  SubtractTowersCS.C \
  SubtractTowersCS_Dict.C \
  RetowerCEMC.C \
  RetowerCEMC_Dict.C \
  CopyAndSubtractJets.C \
  CopyAndSubtractJets_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests
BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_jetbackground_io \
  testexternals_jetbackground


testexternals_jetbackground_io_SOURCES = testexternals.C
testexternals_jetbackground_io_LDADD = libjetbackground_io.la

testexternals_jetbackground_SOURCES = testexternals.C
testexternals_jetbackground_LDADD = libjetbackground.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
