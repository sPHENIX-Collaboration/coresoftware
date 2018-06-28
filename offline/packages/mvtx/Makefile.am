AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmvtx.la \
  libmvtx_io.la


AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

# sources for mvtx library
libmvtx_la_SOURCES = \
  MvtxClusterizer.cc \
  MvtxClusterizer_Dict.cc \
  MvtxDefs.cc \
  MvtxDefs_Dict.cc 

libmvtx_la_LIBADD = \
  -lfun4all \
  -ltrack_util \
  -ltrack_io

# sources for io library
libmvtx_io_la_SOURCES = \
  MvtxHit.cc \
  MvtxHit_Dict.cc

libmvtx_io_la_LIBADD = \
  libmvtx.la \
  -lfun4all \
  -ltrack_util \
  -ltrack_io



# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  MvtxClusterizer.h \
  MvtxDefs.h \
  MvtxHit.h



################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_mvtx_io \
  testexternals_mvtx

testexternals_mvtx_io_SOURCES = testexternals.cc
testexternals_mvtx_io_LDADD = libmvtx_io.la

testexternals_mvtx_SOURCES = testexternals.cc
testexternals_mvtx_LDADD = libmvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
