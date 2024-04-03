AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libcentrality_io.la \
   libcentrality.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libcentrality_io_la_LIBADD = \
  -lphool

libcentrality_la_LIBADD = \
  libcentrality_io.la \
  -lfun4all \
  -lcalo_io \
  -lmbd_io \
  -lphool \
  -lcdbobjects \
  -lffamodules

pkginclude_HEADERS = \
  CentralityInfo.h \
  CentralityInfov1.h \
  CentralityInfov2.h \
  CentralityReco.h


ROOTDICTS = \
  CentralityInfo_Dict.cc \
  CentralityInfov1_Dict.cc \
  CentralityInfov2_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CentralityInfo_Dict_rdict.pcm \
  CentralityInfov1_Dict_rdict.pcm \
  CentralityInfov2_Dict_rdict.pcm

libcentrality_io_la_SOURCES = \
  $(ROOTDICTS) \
  CentralityInfov1.cc \
  CentralityInfov2.cc

libcentrality_la_SOURCES = \
  CentralityReco.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_centrality_io \
  testexternals_centrality

testexternals_centrality_SOURCES = testexternals.cc
testexternals_centrality_LDADD = libcentrality.la

testexternals_centrality_io_SOURCES = testexternals.cc
testexternals_centrality_io_LDADD = libcentrality_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
