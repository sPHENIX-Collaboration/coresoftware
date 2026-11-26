AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libkfparticle_sphenix_io.la \
  libkfparticle_sphenix.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DHomogeneousField

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

pkginclude_HEADERS = \
  KFParticle_caloInfo.h \
  KFParticle_Config.h \
  KFParticle_Container.h \
  KFParticle_detectorInfo.h \
  KFParticle_DST.h \
  KFParticle_eventData.h \
  KFParticle_eventReconstruction.h \
  KFParticle_nTuple.h \
  KFParticle_sPHENIX.h \
  KFParticle_Tools.h \
  KFParticle_triggerInfo.h \
  KFParticle_truthInfo.h \
  KFParticle_util.h

ROOTDICTS = \
  KFParticle_Container_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libkfparticle_sphenix_io_la_SOURCES = \
  $(ROOTDICTS) \
  KFParticle_Container.cc

libkfparticle_sphenix_la_SOURCES = \
  KFParticle_caloInfo.cc \
  KFParticle_detectorInfo.cc \
  KFParticle_DST.cc \
  KFParticle_eventData.cc \
  KFParticle_eventReconstruction.cc \
  KFParticle_nTuple.cc \
  KFParticle_sPHENIX.cc \
  KFParticle_Tools.cc \
  KFParticle_triggerInfo.cc \
  KFParticle_truthInfo.cc \
  KFParticle_util.cc

libkfparticle_sphenix_io_la_LIBADD = \
  -lcentrality \
  -lKFParticle \
  -lphool

libkfparticle_sphenix_la_LIBADD = \
  libkfparticle_sphenix_io.la \
  -lcalo_io \
  -lcalotrigger \
  -lCLHEP \
  -lffarawobjects \
  -lfun4all \
  -lg4eval \
  -lphhepmc \
  -lqautils

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
