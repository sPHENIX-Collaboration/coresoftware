AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libtracker.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libtracker_la_LIBADD = \
  -lCLHEP \
  -lphool \
  -lg4detectors

pkginclude_HEADERS = \
  TrkrCluster.h \
  TrkrClusterv1.h \
  TrkrClusterContainer.h \
  TrackerDefs.h \
  TrkrHitSet.h \
  TrkrHitSetv1.h \
  TrkrHitSetContainer.h

libtracker_la_SOURCES = \
  TrkrCluster.C \
  TrkrCluster_Dict.C \
  TrkrClusterv1.C \
  TrkrClusterv1_Dict.C \
  TrkrClusterContainer.C \
  TrkrClusterContainer_Dict.C \
  TrackerDefs.C \
  TrackerDefs_Dict.C \
  TrkrHitSet.C \
  TrkrHitSet_Dict.C \
  TrkrHitSetv1.C \
  TrkrHitSetv1_Dict.C \
  TrkrHitSetContainer.C \
  TrkrHitSetContainer_Dict.C 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.C
testexternals_tracker_LDADD = libtracker.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
