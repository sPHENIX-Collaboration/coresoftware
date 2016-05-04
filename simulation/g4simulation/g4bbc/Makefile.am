AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libg4bbc_io.la \
   libg4bbc.la

AM_CXXFLAGS = \
 -Wall -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4bbc_io_la_LDFLAGS = \
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
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

%_Dict.cpp: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4bbc_io.la \
  libg4bbc.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
