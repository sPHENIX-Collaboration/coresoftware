AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libcalotrigger_io.la \
   libcalotrigger.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalotrigger_io_la_LIBADD = \
  -lphool

libcalotrigger_la_LIBADD = \
  libcalotrigger_io.la \
  -lcalo_io \
  -lSubsysReco

pkginclude_HEADERS = \
  CaloTriggerSim.h \
  CaloTriggerInfo.h \
  CaloTriggerInfov1.h

ROOTDICTS = \
  CaloTriggerInfo_Dict.cc \
  CaloTriggerInfov1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CaloTriggerInfo_Dict_rdict.pcm \
  CaloTriggerInfov1_Dict_rdict.pcm

libcalotrigger_io_la_SOURCES = \
  $(ROOTDICTS) \
  CaloTriggerInfov1.cc

libcalotrigger_la_SOURCES = \
  CaloTriggerSim.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_calotrigger_io \
  testexternals_calotrigger


testexternals_calotrigger_io_SOURCES = testexternals.cc
testexternals_calotrigger_io_LDADD = libcalotrigger_io.la

testexternals_calotrigger_SOURCES = testexternals.cc
testexternals_calotrigger_LDADD = libcalotrigger.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
