AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_util.la \
  libtrack_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

# sources for util library
libtrack_util_la_SOURCES = \
  TrkrDefs.cc \
  TrkrDefs_Dict.cc

# sources for io library
libtrack_io_la_SOURCES = \
  TrkrCluster_Dict.cc \
  TrkrClusterv1.cc \
  TrkrClusterv1_Dict.cc \
  TrkrClusterContainer.cc \
  TrkrClusterContainer_Dict.cc \
  TrkrClusterHitAssoc.cc \
  TrkrClusterHitAssoc_Dict.cc \
  TrkrHit_Dict.cc \
  TrkrHitTruthAssoc.cc \
  TrkrHitTruthAssoc_Dict.cc \
  TrkrHitSet.cc \
  TrkrHitSet_Dict.cc \
  TrkrHitSetContainer.cc \
  TrkrHitSetContainer_Dict.cc \
  TrkrTrack_Dict.cc \
  TrkrTrackContainer_Dict.cc


libtrack_util_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtrack_util_la_LIBADD = \
  -lCLHEP \
  -lphool


libtrack_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtrack_io_la_LIBADD = \
  -lCLHEP \
  libtrack_util.la


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  TrkrCluster.h \
  TrkrClusterv1.h \
  TrkrClusterContainer.h \
  TrkrClusterHitAssoc.h \
  TrkrDefs.h \
  TrkrHit.h \
  TrkrHitTruthAssoc.h \
  TrkrHitSet.h \
  TrkrHitSetContainer.h \
  TrkrTrack.h \
  TrkrTrackContainer.h


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_track_util \
  testexternals_track_io

testexternals_track_util_SOURCES = testexternals.cc
testexternals_track_util_LDADD = libtrack_util.la

testexternals_track_io_SOURCES = testexternals.cc
testexternals_track_io_LDADD = libtrack_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
