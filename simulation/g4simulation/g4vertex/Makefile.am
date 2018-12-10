AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4vertex_io.la \
   libg4vertex.la

AM_CXXFLAGS = \
 -Werror -msse2

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
  -lg4hough_io \
  libg4vertex_io.la

pkginclude_HEADERS = \
  GlobalVertex.h \
  GlobalVertex_v1.h \
  GlobalVertexFastSimReco.h \
  GlobalVertexMap.h \
  GlobalVertexMap_v1.h \
  GlobalVertexReco.h

libg4vertex_io_la_SOURCES = \
  GlobalVertex.C \
  GlobalVertex_Dict.C \
  GlobalVertex_v1.C \
  GlobalVertex_v1_Dict.C \
  GlobalVertexMap.C \
  GlobalVertexMap_Dict.C \
  GlobalVertexMap_v1.C \
  GlobalVertexMap_v1_Dict.C

libg4vertex_la_SOURCES = \
  GlobalVertexFastSimReco.C \
  GlobalVertexFastSimReco_Dict.C \
  GlobalVertexReco.C \
  GlobalVertexReco_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4vertex_io \
  testexternals_g4vertex

testexternals_g4vertex_io_SOURCES = testexternals.C
testexternals_g4vertex_io_LDADD = libg4vertex_io.la

testexternals_g4vertex_SOURCES = testexternals.C
testexternals_g4vertex_LDADD = libg4vertex.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
