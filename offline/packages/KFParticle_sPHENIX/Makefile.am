AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libkfparticle_sphenix.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -DHomogeneousField
  

pkginclude_HEADERS = \
  KFParticle_Container.h \
  KFParticle_DST.h \
  KFParticle_truthAndDetTools.h \
  KFParticle_nTuple.h \
  KFParticle_Tools.h \
  KFParticle_particleList.h \
  KFParticle_MVA.h \
  KFParticle_eventReconstruction.h \
  KFParticle_sPHENIX.h

if ! MAKEROOT6
  ROOT5_DICTS = \
    KFParticle_Container_DICT.cc \
    KFParticle_DST_DICT.cc \
    KFParticle_truthAndDetTools_Dict.cc \
    KFParticle_nTuple_Dict.cc \
    KFParticle_Tools_Dict.cc \
    KFParticle_particleList_Dict.cc \
    KFParticle_MVA_Dict.cc \
    KFParticle_eventReconstruction_Dict.cc \
    KFParticle_sPHENIX_Dict.cc
endif

libkfparticle_sphenix_la_SOURCES = \
  $(ROOT5_DICTS) \
  KFParticle_Container.cc \
  KFParticle_DST.cc \
  KFParticle_truthAndDetTools.cc \
  KFParticle_nTuple.cc \
  KFParticle_Tools.cc \
  KFParticle_particleList.cc \
  KFParticle_MVA.cc \
  KFParticle_eventReconstruction.cc \
  KFParticle_sPHENIX.cc

libkfparticle_sphenix_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lfun4all \
  -lg4detectors_io \
  -lphg4hit \
  -lg4dst \
  -lg4eval \
  -lKFParticle \
  -L$(ROOTSYS)/lib `root-config --libs` \
  -lTMVA

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libkfparticle_sphenix.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
