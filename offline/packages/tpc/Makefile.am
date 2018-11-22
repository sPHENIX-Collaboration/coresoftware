AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libtpc.la \
  libtpc_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

# sources for tpc library
libtpc_la_SOURCES = \
  TpcDefs.cc \
  TpcDefs_Dict.cc \
  TPCClusterizer.cc \
  TPCClusterizer_Dict.cc 

libtpc_la_LIBADD = \
  -lfun4all \
  -lg4detectors \
  -lSpectrum \
  -ltrack_io \
  -lg4mvtx

# sources for io library
libtpc_io_la_SOURCES = \
  TpcHit.cc \
  TpcHit_Dict.cc

libtpc_io_la_LIBADD = \
  libtpc.la \
  -lfun4all \
  -ltrack_util \
  -ltrack_io \
  -lg4hough_io \
  -lg4detectors_io \
  -ltrackbase_historic_io


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  TpcDefs.h \
  TpcHit.h \
  TPCClusterizer.h

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_tpc_io \
  testexternals_tpc

testexternals_tpc_io_SOURCES = testexternals.cc
testexternals_tpc_io_LDADD = libtpc_io.la

testexternals_tpc_SOURCES = testexternals.cc
testexternals_tpc_LDADD = libtpc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
