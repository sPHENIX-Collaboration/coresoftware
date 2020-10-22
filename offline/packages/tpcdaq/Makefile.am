##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libtpcdaq.la \
   libtpcdaq_io.la

libtpcdaq_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtpcdaq_la_LIBADD = \
  -lfun4all \
  -lg4dst \
  -lphool

libtpcdaq_io_la_LIBADD = \
  -lphool

pkginclude_HEADERS = \
  TPCDaqDefs.h

ROOTDICTS = \
  TPCFEETestRecov1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TPCFEETestRecov1_Dict_rdict.pcm

libtpcdaq_la_SOURCES = \
  TPCIntegratedCharge.cc \
  TPCDataStreamEmulator.cc \
  TPCDaqDefs.cc \
  TPCDaqDefs_Dict.cc

libtpcdaq_io_la_SOURCES = \
  $(ROOTDICTS) \
  TPCFEETestRecov1.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals \
  testexternals_io

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libtpcdaq.la

testexternals_io_SOURCES = \
  testexternals.cc

testexternals_io_LDADD = \
  libtpcdaq_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
