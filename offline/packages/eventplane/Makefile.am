AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libeventplane_io.la \
  libeventplane.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  EpFinder.h \
  EpInfo.h \
  EpInfov1.h \
  EpFinderReco.h

ROOTDICTS = \
  EpInfo_Dict.cc \
  EpInfov1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  EpInfo_Dict_rdict.pcm \
  EpInfov1_Dict_rdict.pcm

libeventplane_io_la_SOURCES = \
  $(ROOTDICTS) \
  EpInfov1.cc

libeventplane_la_SOURCES = \
  EpFinder.cc \
  EpFinderReco.cc

libeventplane_io_la_LIBADD = \
  -lphool

libeventplane_la_LIBADD = \
  libeventplane_io.la \
  -lphool \
  -lSubsysReco \
  -lcentrality_io \
  -lphg4hit \
  -lgsl \
  -lgslcblas \
  -lcalo_io

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

#################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libeventplane_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libeventplane.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
