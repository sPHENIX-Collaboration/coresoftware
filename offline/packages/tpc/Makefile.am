##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libtpc_io.la \
  libtpc.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include/eigen3

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64
pkginclude_HEADERS = \
  TpcDefs.h \
  TpcHit.h \
  TpcClusterizer.h \
  TpcSpaceChargeCorrection.h \
  TpcClusterCleaner.h

ROOTDICTS = \
  TpcHit_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
    TpcHit_Dict_rdict.pcm

# sources for tpc library
libtpc_la_SOURCES = \
  TpcClusterizer.cc \
  TpcSpaceChargeCorrection.cc \
  TpcClusterCleaner.cc

libtpc_la_LIBADD = \
  libtpc_io.la \
  -lfun4all \
  -lSpectrum \
  -lpthread

# sources for io library
libtpc_io_la_SOURCES = \
  $(ROOTDICTS) \
  TpcDefs.cc \
  TpcHit.cc

libtpc_io_la_LIBADD = \
  -ltrack_io \
  -lg4detectors_io \
  -ltrackbase_historic_io


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

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
	rm -f *Dict* *.pcm $(BUILT_SOURCES)
