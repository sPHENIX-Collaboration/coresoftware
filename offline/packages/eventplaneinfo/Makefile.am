AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libeventplaneinfo_io.la \
   libeventplaneinfo.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libeventplaneinfo_io_la_LIBADD = \
  -lphool

libeventplaneinfo_la_LIBADD = \
  libeventplaneinfo_io.la \
  -lepd_io \
  -lmbd_io \
  -lfun4all \
  -ltrackbase_historic_io

pkginclude_HEADERS = \
  Eventplaneinfo.h \
  Eventplaneinfov1.h \
  EventplaneinfoMap.h \
  EventplaneinfoMapv1.h \
  EventPlaneReco.h

ROOTDICTS = \
  Eventplaneinfo_Dict.cc \
  Eventplaneinfov1_Dict.cc \
  EventplaneinfoMap_Dict.cc \
  EventplaneinfoMapv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  Eventplaneinfo_Dict_rdict.pcm \
  Eventplaneinfov1_Dict_rdict.pcm \
  EventplaneinfoMap_Dict_rdict.pcm \
  EventplaneinfoMapv1_Dict_rdict.pcm

libeventplaneinfo_io_la_SOURCES = \
  $(ROOTDICTS) \
  Eventplaneinfo.cc \
  Eventplaneinfov1.cc \
  EventplaneinfoMap.cc \
  EventplaneinfoMapv1.cc

libeventplaneinfo_la_SOURCES = \
  EventPlaneReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_eventplaneinfo_io \
  testexternals_eventplaneinfo

BUILT_SOURCES = \
  testexternals.cc

testexternals_eventplaneinfo_io_SOURCES = testexternals.cc
testexternals_eventplaneinfo_io_LDADD = libeventplaneinfo_io.la

testexternals_eventplaneinfo_SOURCES = testexternals.cc
testexternals_eventplaneinfo_LDADD = libeventplaneinfo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
