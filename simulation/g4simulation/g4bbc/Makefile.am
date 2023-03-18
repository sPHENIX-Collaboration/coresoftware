AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libg4bbc_io.la \
   libg4bbc.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4bbc_io_la_LIBADD = \
  -lphool

libg4bbc_la_LIBADD = \
  libg4bbc_io.la \
  -lg4detectors \
  -lfun4all \
  -lSubsysReco \
  -lbbc

pkginclude_HEADERS = \
  BbcVertex.h \
  BbcVertexv1.h \
  BbcVertexMap.h \
  BbcVertexMapv1.h \
  BbcSimReco.h \
  BbcVertexFastSimReco.h

ROOTDICTS = \
  BbcVertex_Dict.cc \
  BbcVertexv1_Dict.cc \
  BbcVertexMap_Dict.cc \
  BbcVertexMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  BbcVertex_Dict_rdict.pcm \
  BbcVertexv1_Dict_rdict.pcm \
  BbcVertexMap_Dict_rdict.pcm \
  BbcVertexMapv1_Dict_rdict.pcm

libg4bbc_io_la_SOURCES = \
  $(ROOTDICTS) \
  BbcVertexv1.cc \
  BbcVertexMap.cc \
  BbcVertexMapv1.cc

libg4bbc_la_SOURCES = \
  BbcSimReco.cc \
  BbcVertexFastSimReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternalsg4bbc_io \
  testexternalsg4bbc

BUILT_SOURCES = testexternals.cc

testexternalsg4bbc_io_SOURCES = testexternals.cc
testexternalsg4bbc_io_LDADD = libg4bbc_io.la

testexternalsg4bbc_SOURCES = testexternals.cc
testexternalsg4bbc_LDADD = libg4bbc_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
