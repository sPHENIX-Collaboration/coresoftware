##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib \
  -L$(OPT_SPHENIX)/lib

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
  CylinderGeomInttHelper.h \
  InttArborist.h \
  InttBadChannelMap.h \
  InttBCOMap.h \
  InttClusterizer.h \
  InttCombinedRawDataDecoder.h \
  InttDacMap.h \
  InttFeeMap.h \
  InttFelixMap.h \
  InttMap.h \
  InttMapping.h \
  InttOdbcQuery.h \
  InttSurveyMap.h \
  InttXYVertexFinder.h \
  InttZVertexFinder.h \
  InttVertexMap.h \
  InttVertexMapv1.h \
  InttVertex.h \
  InttVertexv1.h

noinst_HEADERS = \
  INTTZvtx.h \
  INTTXYvtx.h \
  InttVertexUtil.h

ROOTDICTS = \
  CylinderGeomIntt_Dict.cc \
  InttVertexMap_Dict.cc \
  InttVertexMapv1_Dict.cc \
  InttVertex_Dict.cc \
  InttVertexv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CylinderGeomIntt_Dict_rdict.pcm \
  InttVertexMap_Dict_rdict.pcm \
  InttVertexMapv1_Dict_rdict.pcm \
  InttVertex_Dict_rdict.pcm \
  InttVertexv1_Dict_rdict.pcm

# sources for intt library
libintt_la_SOURCES = \
  CylinderGeomInttHelper.cc \
  InttArborist.cc \
  InttBCOMap.cc \
  InttBadChannelMap.cc \
  InttClusterizer.cc \
  InttCombinedRawDataDecoder.cc \
  InttDacMap.cc \
  InttFeeMap.cc \
  InttFelixMap.cc \
  InttMapping.cc \
  InttOdbcQuery.cc \
  InttMap.cc \
  InttSurveyMap.cc \
  InttVertexUtil.cc \
  InttXYVertexFinder.cc \
  InttZVertexFinder.cc \
  INTTXYvtx.cc \
  INTTZvtx.cc

libintt_la_LIBADD = \
  libintt_io.la \
  -lActsCore \
  -lcdbobjects \
  -lCLHEP \
  -lffarawobjects \
  -lodbc++ \
  -lphg4hit \
  -lSubsysReco \
  -ltrack

# sources for io library
libintt_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeomIntt.cc \
  InttVertexMap.cc \
  InttVertexMapv1.cc \
  InttVertexv1.cc

libintt_io_la_LIBADD = \
  -lffamodules \
  -lg4detectors_io \
  -lphool \
  -lglobalvertex_io

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
