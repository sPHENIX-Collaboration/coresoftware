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
  TrkrDefUtil.C \
  TrkrDefUtil_Dict.C

# sources for io library
libtrack_io_la_SOURCES = \
  TrkrCluster.C \
  TrkrCluster_Dict.C \
  TrkrClusterv1.C \
  TrkrClusterv1_Dict.C \
  TrkrClusterContainer.C \
  TrkrClusterContainer_Dict.C \
  TrkrHitSet.C \
  TrkrHitSet_Dict.C \
  TrkrHitSetv1.C \
  TrkrHitSetv1_Dict.C \
  TrkrHitSetContainer.C \
  TrkrHitSetContainer_Dict.C 


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


#libtracker_la_LIBADD = \
#  -lCLHEP \
#  -lphool \
#  -lg4detectors


# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  TrkrCluster.h \
  TrkrClusterv1.h \
  TrkrClusterContainer.h \
  TrkrDefUtil.h \
  TrkrHitSet.h \
  TrkrHitSetv1.h \
  TrkrHitSetContainer.h


################################################
# linking tests

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_track_util \
  testexternals_track_io

testexternals_track_util_SOURCES = testexternals.C
testexternals_track_util_LDADD = libtrack_util.la

testexternals_track_io_SOURCES = testexternals.C
testexternals_track_io_LDADD = libtrack_io.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
