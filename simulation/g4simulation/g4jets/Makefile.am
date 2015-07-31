AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4 \
  -I`root-config --incdir`

lib_LTLIBRARIES = libg4jets.la

AM_CXXFLAGS = \
 -Wall -Werror -msse2

libg4jets_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lfun4all \
  -lphg4hit \
  -lg4hough_io \
  -lcemc_io \
  -lCGAL \
  -lfastjet 

pkginclude_HEADERS = \
  Jet.h \
  JetV1.h \
  JetMap.h \
  JetMapV1.h \
  JetInput.h \
  JetAlgo.h \
  JetReco.h \
  TruthJetInput.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4jets_la_SOURCES = \
  Jet.C \
  Jet_Dict.C \
  JetV1.C \
  JetV1_Dict.C \
  JetMap.C \
  JetMap_Dict.C \
  JetMapV1.C \
  JetMapV1_Dict.C \
  JetReco.C \
  JetReco_Dict.C \
  TruthJetInput.C \
  TruthJetInput_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

%_Dict.cpp: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4jets.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
