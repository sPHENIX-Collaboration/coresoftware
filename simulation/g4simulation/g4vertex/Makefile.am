AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4vertex_io.la \
   libg4vertex.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4vertex_io_la_LIBADD = \
  -lphool \
  -ltrackbase_historic_io

libg4vertex_la_LIBADD = \
  -lg4detectors \
  -lfun4all \
  -lg4bbc_io \
  -ltrackbase_historic_io \
  libg4vertex_io.la

pkginclude_HEADERS = \
  GlobalVertex.h \
  GlobalVertexv1.h \
  GlobalVertexFastSimReco.h \
  GlobalVertexMap.h \
  GlobalVertexMapv1.h \
  GlobalVertexReco.h

ROOTDICTS = \
  GlobalVertex_Dict.cc \
  GlobalVertexv1_Dict.cc \
  GlobalVertexMap_Dict.cc \
  GlobalVertexMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  GlobalVertex_Dict_rdict.pcm \
  GlobalVertexv1_Dict_rdict.pcm \
  GlobalVertexMap_Dict_rdict.pcm \
  GlobalVertexMapv1_Dict_rdict.pcm

libg4vertex_io_la_SOURCES = \
  $(ROOTDICTS) \
  GlobalVertex.cc \
  GlobalVertexv1.cc \
  GlobalVertexMap.cc \
  GlobalVertexMapv1.cc

libg4vertex_la_SOURCES = \
  GlobalVertexFastSimReco.cc \
  GlobalVertexReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4vertex_io \
  testexternals_g4vertex

BUILT_SOURCES = \
  testexternals.cc

testexternals_g4vertex_io_SOURCES = testexternals.cc
testexternals_g4vertex_io_LDADD = libg4vertex_io.la

testexternals_g4vertex_SOURCES = testexternals.cc
testexternals_g4vertex_LDADD = libg4vertex.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
