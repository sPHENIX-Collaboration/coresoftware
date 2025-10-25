AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  CaloVtxReco.h

ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

lib_LTLIBRARIES = \
  libcalovtxreco_io.la \
  libcalovtxreco.la

libcalovtxreco_io_la_SOURCES = \
  $(ROOTDICTS) 

libcalovtxreco_la_SOURCES = \
  CaloVtxReco.cc

libcalovtxreco_io_la_LIBADD = \
  -lphool -lSubsysReco -lmbd_io -ljetbase

libcalovtxreco_la_LIBADD = \
  libcalovtxreco_io.la \
  -lcalo_io \
  -lSubsysReco -lffaobjects \
  -ljetbase -lglobalvertex
#-lg4jets_io
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libcalovtxreco_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcalovtxreco.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
