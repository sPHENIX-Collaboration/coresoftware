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
  KFParticle_truthAndDetTools.h \
  KFParticle_nTuple.h \
  KFParticle_Tools.h \
  KFParticle_particleList.h \
  KFParticle_MVA.h \
  KFParticle_sPHENIX.h

if ! MAKEROOT6
  ROOT5_DICTS = \
    KFParticle_truthAndDetTools_Dict.cc \
    KFParticle_nTuple_Dict.cc \
    KFParticle_Tools_Dict.cc \
    KFParticle_particleList_Dict.cc \
    KFParticle_MVA_Dict.cc \
    KFParticle_sPHENIX_Dict.cc
endif

libkfparticle_sphenix_la_SOURCES = \
  $(ROOT5_DICTS) \
  KFParticle_truthAndDetTools.cxx \
  KFParticle_nTuple.cxx \
  KFParticle_Tools.cxx \
  KFParticle_particleList.cxx \
  KFParticle_MVA.cxx \
  KFParticle_sPHENIX.cxx

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
