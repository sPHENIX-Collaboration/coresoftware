AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

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
  MbdVertex.h \
  MbdVertexv1.h \
  MbdVertexv2.h \
  MbdVertexMap.h \
  MbdVertexMapv1.h \
  SvtxVertex.h \
  SvtxVertex_v1.h \
  SvtxVertex_v2.h \
  SvtxVertexMap.h \
  SvtxVertexMap_v1.h \
  Vertex.h

ROOTDICTS = \
  GlobalVertex_Dict.cc \
  GlobalVertexv1_Dict.cc \
  GlobalVertexv2_Dict.cc \
  GlobalVertexMap_Dict.cc \
  GlobalVertexMapv1_Dict.cc \
  MbdVertex_Dict.cc \
  MbdVertexv1_Dict.cc \
  MbdVertexv2_Dict.cc \
  MbdVertexMap_Dict.cc \
  MbdVertexMapv1_Dict.cc \
  SvtxVertex_Dict.cc \
  SvtxVertex_v1_Dict.cc \
  SvtxVertex_v2_Dict.cc \
  SvtxVertexMap_Dict.cc \
  SvtxVertexMap_v1_Dict.cc \
  Vertex_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libglobalvertex_io_la_SOURCES = \
  $(ROOTDICTS) \
  GlobalVertex.cc \
  GlobalVertexv1.cc \
  GlobalVertexv2.cc \
  GlobalVertexMap.cc \
  GlobalVertexMapv1.cc \
  MbdVertexv1.cc \
  MbdVertexv2.cc \
  MbdVertexMap.cc \
  MbdVertexMapv1.cc \
  SvtxVertex.cc \
  SvtxVertex_v1.cc \
  SvtxVertex_v2.cc \
  SvtxVertexMap.cc \
  SvtxVertexMap_v1.cc \
  Vertex.cc

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
