AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmvtx_io.la \
  libmvtx.la 


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
  MVTXClusterizer.cc \
  MVTXClusterizer_Dict.cc \
  MvtxClusterizer.cc \
  MvtxClusterizer_Dict.cc \
  MvtxDefs.cc \
  MvtxDefs_Dict.cc

libmvtx_la_LIBADD = \
  libmvtx_io.la \
  -lfun4all \
  -ltrack_util \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lg4detectors \
  -lg4mvtx 

# sources for io library
libmvtx_io_la_SOURCES = \
  MvtxHit.cc \
  MvtxHit_Dict.cc

libmvtx_io_la_LIBADD = \
  -lphool \
  -ltrack_io

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  MVTXClusterizer.h \
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
