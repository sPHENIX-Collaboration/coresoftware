AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libglobalvertex_io.la \
   libglobalvertex.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libglobalvertex_io_la_LIBADD = \
  -lphool \
  -ltrackbase_historic_io

libglobalvertex_la_LIBADD = \
  libglobalvertex_io.la \
  -ltrackbase_historic_io \
  -lg4detectors \
  -lfun4all

pkginclude_HEADERS = \
  GlobalVertex.h \
  GlobalVertexv1.h \
  GlobalVertexv2.h \
  GlobalVertexMap.h \
  GlobalVertexMapv1.h \
  GlobalVertexReco.h \
  SvtxVertex.h \
  SvtxVertex_v1.h \
  SvtxVertex_v2.h \
  SvtxVertexMap.h \
  SvtxVertexMap_v1.h \
  MbdVertex.h \
  MbdVertexv1.h \
  MbdVertexv2.h \
  MbdVertexMap.h \
  MbdVertexMapv1.h \
  Vertex.h

ROOTDICTS = \
  Vertex_Dict.cc \
  GlobalVertex_Dict.cc \
  GlobalVertexv1_Dict.cc \
  GlobalVertexv2_Dict.cc \
  GlobalVertexMap_Dict.cc \
  GlobalVertexMapv1_Dict.cc \
  SvtxVertex_Dict.cc \
  SvtxVertex_v1_Dict.cc \
  SvtxVertex_v2_Dict.cc \
  SvtxVertexMap_Dict.cc \
  SvtxVertexMap_v1_Dict.cc \
  MbdVertex_Dict.cc \
  MbdVertexv1_Dict.cc \
  MbdVertexv2_Dict.cc \
  MbdVertexMap_Dict.cc \
  MbdVertexMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  Vertex_Dict_rdict.pcm \
  GlobalVertex_Dict_rdict.pcm \
  GlobalVertexv1_Dict_rdict.pcm \
  GlobalVertexv2_Dict_rdict.pcm \
  GlobalVertexMap_Dict_rdict.pcm \
  GlobalVertexMapv1_Dict_rdict.pcm \
  SvtxVertex_Dict_rdict.pcm \
  SvtxVertex_v1_Dict_rdict.pcm \
  SvtxVertex_v2_Dict_rdict.pcm \
  SvtxVertexMap_Dict_rdict.pcm \
  SvtxVertexMap_v1_Dict_rdict.pcm \
  MbdVertex_Dict_rdict.pcm \
  MbdVertexv1_Dict_rdict.pcm \
  MbdVertexv2_Dict_rdict.pcm \
  MbdVertexMap_Dict_rdict.pcm \
  MbdVertexMapv1_Dict_rdict.pcm

libglobalvertex_io_la_SOURCES = \
  $(ROOTDICTS) \
  Vertex.cc \
  GlobalVertex.cc \
  GlobalVertexv1.cc \
  GlobalVertexv2.cc \
  GlobalVertexMap.cc \
  GlobalVertexMapv1.cc \
  SvtxVertex.cc \
  SvtxVertex_v1.cc \
  SvtxVertex_v2.cc \
  SvtxVertexMap.cc \
  SvtxVertexMap_v1.cc \
  MbdVertexv1.cc \
  MbdVertexv2.cc \
  MbdVertexMap.cc \
  MbdVertexMapv1.cc

libglobalvertex_la_SOURCES = \
  GlobalVertexReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_globalvertex_io \
  testexternals_globalvertex

BUILT_SOURCES = \
  testexternals.cc

testexternals_globalvertex_io_SOURCES = testexternals.cc
testexternals_globalvertex_io_LDADD = libglobalvertex_io.la

testexternals_globalvertex_SOURCES = testexternals.cc
testexternals_globalvertex_LDADD = libglobalvertex.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
