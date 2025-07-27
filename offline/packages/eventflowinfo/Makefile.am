AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libeventflowinfo_io.la \
   libeventflowinfo.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libeventflowinfo_io_la_LIBADD = \
  -lphool \
  -lcentrality_io 

libeventflowinfo_la_LIBADD = \
  libeventflowinfo_io.la \
  -lepd_io \
  -lmbd_io \
  -lffamodules \
  -lcalotrigger_io \
  -lcalo_io \
  -lcentrality_io \
  -ljetbase \
  -lfun4all \
  -lffarawobjects \
  -lglobalvertex_io

pkginclude_HEADERS = \
  EventFlowInfo.h \
  EventFlowInfov1.h \
  EventFlowInfoMap.h \
  EventFlowInfoMapv1.h \
  EventFlowReco.h

ROOTDICTS = \
  EventFlowInfo_Dict.cc \
  EventFlowInfov1_Dict.cc \
  EventFlowInfoMap_Dict.cc \
  EventFlowInfoMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  EventFlowInfov1_Dict_rdict.pcm \
  EventFlowInfoMapv1_Dict_rdict.pcm

libeventflowinfo_io_la_SOURCES = \
  $(ROOTDICTS) \
  EventFlowInfov1.cc \
  EventFlowInfoMapv1.cc

libeventflowinfo_la_SOURCES = \
  EventFlowReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_eventflowinfo_io \
  testexternals_eventflowinfo

BUILT_SOURCES = \
  testexternals.cc

testexternals_eventflowinfo_io_SOURCES = testexternals.cc
testexternals_eventflowinfo_io_LDADD = libeventflowinfo_io.la

testexternals_eventflowinfo_SOURCES = testexternals.cc
testexternals_eventflowinfo_LDADD = libeventflowinfo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
