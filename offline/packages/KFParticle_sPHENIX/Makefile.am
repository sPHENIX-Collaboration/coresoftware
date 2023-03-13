AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libkfparticle_sphenix_io.la \
    libkfparticle_sphenix.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField


pkginclude_HEADERS = \
  KFParticle_Container.h \
  KFParticle_DST.h \
  KFParticle_truthAndDetTools.h \
  KFParticle_nTuple.h \
  KFParticle_Tools.h \
  KFParticle_MVA.h \
  KFParticle_eventReconstruction.h \
  KFParticle_sPHENIX.h

ROOTDICTS = \
  KFParticle_Container_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  KFParticle_Container_Dict_rdict.pcm

libkfparticle_sphenix_io_la_SOURCES = \
  $(ROOTDICTS) \
  KFParticle_Container.cc

libkfparticle_sphenix_la_SOURCES = \
  KFParticle_DST.cc \
  KFParticle_truthAndDetTools.cc \
  KFParticle_nTuple.cc \
  KFParticle_Tools.cc \
  KFParticle_MVA.cc \
  KFParticle_eventReconstruction.cc \
  KFParticle_sPHENIX.cc

libkfparticle_sphenix_io_la_LIBADD = \
  -lKFParticle \
  -lphool

libkfparticle_sphenix_la_LIBADD = \
  libkfparticle_sphenix_io.la \
  -lcalo_io \
  -lfun4all \
  -lg4eval \
  -lTMVA \
  -lphhepmc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals \
  testexternals_io

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libkfparticle_sphenix.la

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD = libkfparticle_sphenix_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
