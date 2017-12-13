AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libmvtx.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

libmvtx_la_LIBADD = \
  -lCLHEP \
  -lphool \
  -lg4detectors \
  -lg4hough \
  -ltracker

pkginclude_HEADERS = \
  MvtxClusterizer.h \
  MvtxDigitizer.h \
  MvtxThresholds.h

libmvtx_la_SOURCES = \
  MvtxClusterizer.C \
  MvtxClusterizer_Dict.C \
  MvtxDigitizer.C \
  MvtxDigitizer_Dict.C \
  MvtxThresholds.C \
  MvtxThresholds_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_tracker

testexternals_tracker_SOURCES = testexternals.C
testexternals_tracker_LDADD = libmvtx.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm
