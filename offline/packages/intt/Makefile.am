##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

if USE_ONLINE

lib_LTLIBRARIES = \
  libintt.la 

pkginclude_HEADERS = \
  InttFelixMap.h \
  InttMapping.h

libintt_la_SOURCES = \
  InttFelixMap.cc \
  InttMapping.cc

else

lib_LTLIBRARIES = \
  libintt_io.la \
  libintt.la 

pkginclude_HEADERS = \
  CylinderGeomIntt.h \
  InttArborist.h \
  InttBadChannelMap.h \
  InttBadChannelMapv1.h \
  InttBCOMap.h \
  InttClusterizer.h \
  InttCombinedRawDataDecoder.h \
  InttDacMap.h \
  InttFeeMap.h \
  InttFeeMapv1.h \
  InttFelixMap.h \
  InttMap.h \
  InttMapping.h \
  InttSurveyMap.h \
  InttSurveyMapv1.h \
  InttVertexFinder.h \
  InttVertexMap.h \
  InttVertexMapv1.h \
  InttVertex.h \
  InttVertexv1.h

ROOTDICTS = \
  CylinderGeomIntt_Dict.cc \
  InttBadChannelMap_Dict.cc \
  InttBadChannelMapv1_Dict.cc \
  InttFeeMap_Dict.cc \
  InttFeeMapv1_Dict.cc \
  InttMap_Dict.cc \
  InttSurveyMap_Dict.cc \
  InttSurveyMapv1_Dict.cc \
  InttVertexMap_Dict.cc \
  InttVertexMapv1_Dict.cc \
  InttVertex_Dict.cc \
  InttVertexv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CylinderGeomIntt_Dict_rdict.pcm \
  InttBadChannelMap_Dict_rdict.pcm \
  InttBadChannelMapv1_Dict_rdict.pcm \
  InttFeeMap_Dict_rdict.pcm \
  InttFeeMapv1_Dict_rdict.pcm \
  InttMap_Dict_rdict.pcm \
  InttSurveyMap_Dict_rdict.pcm \
  InttSurveyMapv1_Dict_rdict.pcm \
  InttVertexMap_Dict_rdict.pcm \
  InttVertexMapv1_Dict_rdict.pcm \
  InttVertex_Dict_rdict.pcm \
  InttVertexv1_Dict_rdict.pcm

# sources for intt library
libintt_la_SOURCES = \
  InttArborist.cc \
  InttBCOMap.cc \
  InttBadChannelMap.cc \
  InttBadChannelMapv1.cc \
  InttClusterizer.cc \
  InttCombinedRawDataDecoder.cc \
  InttDacMap.cc \
  InttFeeMap.cc \
  InttFeeMapv1.cc \
  InttFelixMap.cc \
  InttMapping.cc \
  InttMap.cc \
  InttSurveyMap.cc \
  InttSurveyMapv1.cc \
  InttVertexFinder.cc

libintt_la_LIBADD = \
  libintt_io.la \
  -lcdbobjects \
  -lCLHEP \
  -lffarawobjects \
  -lphg4hit \
  -lSubsysReco

# sources for io library
libintt_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeomIntt.cc \
  InttBadChannelMap.cc \
  InttBadChannelMapv1.cc \
  InttFeeMap.cc \
  InttFeeMapv1.cc \
  InttMap.cc \
  InttSurveyMap.cc \
  InttSurveyMapv1.cc \
  InttVertexMap.cc \
  InttVertexMapv1.cc \
  InttVertexv1.cc

libintt_io_la_LIBADD = \
  -lcdbobjects \
  -lg4detectors_io \
  -lphool \
  -ltrack_io \
  -ltrackbase_historic_io

endif

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_intt_io \
  testexternals_intt

testexternals_intt_io_SOURCES = testexternals.cc
testexternals_intt_io_LDADD = libintt_io.la

testexternals_intt_SOURCES = testexternals.cc
testexternals_intt_LDADD = libintt.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
