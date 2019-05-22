##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libintt_io.la \
  libintt.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  InttDefs.h \
  InttHit.h \
  InttClusterizer.h \
  CylinderGeomIntt.h

ROOTDICTS = \
  InttHit_Dict.cc \
  CylinderGeomIntt_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  InttHit_Dict_rdict.pcm \
  CylinderGeomIntt_Dict_rdict.pcm
else
  ROOT5_DICTS = \
  InttClusterizer_Dict.cc

endif

# sources for intt library
libintt_la_SOURCES = \
  $(ROOT5_DICTS) \
  InttClusterizer.cc

libintt_la_LIBADD = \
  libintt_io.la \
  -lCLHEP \
  -lphg4hit \
  -lSubsysReco

# sources for io library
libintt_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeomIntt.cc \
  InttDefs.cc \
  InttHit.cc

libintt_io_la_LIBADD = \
  -lCLHEP \
  -lg4detectors_io \
  -lphool \
  -ltrack_io \
  -ltrackbase_historic_io


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_intt_io \
  testexternals_intt

testexternals_intt_io_SOURCES = testexternals.cc
testexternals_intt_io_LDADD = libintt_io.la

testexternals_intt_SOURCES = testexternals.cc
testexternals_intt_LDADD = libintt.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
