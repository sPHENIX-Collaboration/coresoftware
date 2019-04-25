##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmvtx_io.la \
  libmvtx.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  MvtxDefs.h \
  MvtxHit.h \
  MvtxClusterizer.h \
  CylinderGeom_MVTX.h

ROOTDICTS = \
  MvtxHit_Dict.cc \
  CylinderGeom_MVTX_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    MvtxHit_Dict_rdict.pcm \
    CylinderGeom_MVTX_Dict_rdict.pcm
else
  ROOT5_DICTS = \
    MvtxClusterizer_Dict.cc
endif

# sources for mvtx library
libmvtx_la_SOURCES = \
  $(ROOT5_DICTS) \
  MvtxClusterizer.cc

libmvtx_la_LIBADD = \
  libmvtx_io.la \
  -lCLHEP \
  -lphg4hit \
  -lSubsysReco

# sources for io library
libmvtx_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeom_MVTX.cc \
  MvtxDefs.cc \
  MvtxHit.cc

libmvtx_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io \
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
  testexternals_mvtx_io \
  testexternals_mvtx

testexternals_mvtx_io_SOURCES = testexternals.cc
testexternals_mvtx_io_LDADD = libmvtx_io.la

testexternals_mvtx_SOURCES = testexternals.cc
testexternals_mvtx_LDADD = libmvtx.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
