AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I`root-config --incdir` 

lib_LTLIBRARIES = \
   libtracker_io.la \
   libtracker.la

AM_CXXFLAGS = \
 -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtracker_io_la_LIBADD = \
  -lphool

libtracker_la_LDFLAGS = \
  `geant4-config --libs`

libtracker_la_LIBADD = \
  -lg4detectors \
  -lfun4all \
  libtracker_io.la

pkginclude_HEADERS = \
  TrackerCluster.h \
  TrackerDefs.h \
  TrackerHit.h

libtracker_io_la_SOURCES = \
  TrackerCluster.C \
  TrackerCluster_Dict.C \
  TrackerDefs.C \
  TrackerDefs_Dict.C \
  TrackerHit.C \
  TrackerHit_Dict.C

libtracker_la_SOURCES = \

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_tracker_io \
  testexternals_tracker

testexternals_tracker_io_SOURCES = testexternals.C
testexternals_tracker_io_LDADD = libtracker_io.la

testexternals_tracker_SOURCES = testexternals.C
testexternals_tracker_LDADD = libtracker.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* *.pcm
