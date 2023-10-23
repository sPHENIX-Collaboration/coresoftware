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
  -lbbc_io \
  -lg4detectors \
  -lfun4all \
  -lg4bbc_io \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  GlobalVertex.h \
  GlobalVertexv1.h \
  GlobalVertexMap.h \
  GlobalVertexMapv1.h \
  GlobalVertexReco.h \
  Vertex.h

ROOTDICTS = \
  GlobalVertex_Dict.cc \
  GlobalVertexv1_Dict.cc \
  GlobalVertexMap_Dict.cc \
  GlobalVertexMapv1_Dict.cc \
  Vertex_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  GlobalVertex_Dict_rdict.pcm \
  GlobalVertexv1_Dict_rdict.pcm \
  GlobalVertexMap_Dict_rdict.pcm \
  GlobalVertexMapv1_Dict_rdict.pcm \
  Vertex_Dict_rdict.pcm

libglobalvertex_io_la_SOURCES = \
  $(ROOTDICTS) \
  GlobalVertex.cc \
  GlobalVertexv1.cc \
  GlobalVertexMap.cc \
  GlobalVertexMapv1.cc \
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
