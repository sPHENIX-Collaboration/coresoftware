AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmvtx.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libmvtx_la_LIBADD = \
  -lfun4all \
  -ltrack_util \
  -ltrack_io

pkginclude_HEADERS = \
  MvtxClusterizer.h \
  MvtxDefs.h \
  MvtxHit.h

libmvtx_la_SOURCES = \
  MvtxClusterizer.cc \
  MvtxClusterizer_Dict.cc \
  MvtxDefs.cc \
  MvtxDefs_Dict.cc \
  MvtxHit.cc \
  MvtxHit_Dict.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.cc
testexternals_tracker_LDADD = libmvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
