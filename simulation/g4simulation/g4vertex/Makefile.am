AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4vertex.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4vertex_la_LIBADD = \
  -lg4detectors \
  -lfun4all \
  -lg4mbd_io \
  -lmbd_io \
  -ltrackbase_historic_io \
  -lglobalvertex_io

pkginclude_HEADERS = \
  GlobalVertexFastSimReco.h

libg4vertex_la_SOURCES = \
  GlobalVertexFastSimReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4vertex

BUILT_SOURCES = \
  testexternals.cc

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
