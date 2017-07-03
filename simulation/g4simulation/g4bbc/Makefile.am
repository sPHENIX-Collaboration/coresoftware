AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4bbc_io.la \
   libg4bbc.la

AM_CXXFLAGS = \
 -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4bbc_io_la_LIBADD = \
  -lphool

libg4bbc_la_LIBADD = \
  -lg4detectors \
  -lfun4all \
  libg4bbc_io.la

pkginclude_HEADERS = \
  BbcVertex.h \
  BbcVertex_v1.h \
  BbcVertexMap.h \
  BbcVertexMap_v1.h \
  BbcVertexFastSimReco.h

libg4bbc_io_la_SOURCES = \
  BbcVertex.C \
  BbcVertex_Dict.C \
  BbcVertex_v1.C \
  BbcVertex_v1_Dict.C \
  BbcVertexMap.C \
  BbcVertexMap_Dict.C \
  BbcVertexMap_v1.C \
  BbcVertexMap_v1_Dict.C

libg4bbc_la_SOURCES = \
  BbcVertexFastSimReco.C \
  BbcVertexFastSimReco_Dict.C 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternalsg4bbc_io \
  testexternalsg4bbc

testexternalsg4bbc_io_SOURCES = testexternals.C
testexternalsg4bbc_io_LDADD = libg4bbc_io.la

testexternalsg4bbc_SOURCES = testexternals.C
testexternalsg4bbc_LDADD = libg4bbc_io.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
