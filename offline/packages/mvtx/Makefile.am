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
  MVTXClusterizer.h \
  MvtxClusterizer.h \
  MvtxDefs.h \
  MvtxHit.h \
  CylinderGeom_MVTX.h

ROOTDICTS = \
  MvtxHit_Dict.cc
if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    MvtxHit_Dict_rdict.pcm
else
  ROOT5_IO_DICTS = \
    MvtxDefs_Dict.cc

  ROOT5_DICTS = \
    MVTXClusterizer_Dict.cc \
    MvtxClusterizer_Dict.cc \
    CylinderGeom_MVTX_Dict.cc
endif

# sources for io library
libmvtx_io_la_SOURCES = \
  $(ROOTDICTS) \
  $(ROOT5_IO_DICTS) \
  MvtxDefs.cc \
  MvtxHit.cc

libmvtx_io_la_LIBADD = \
  -lphool \
  -ltrack_io

# sources for mvtx library
libmvtx_la_SOURCES = \
  $(ROOT5_DICTS) \
  MVTXClusterizer.cc \
  MvtxClusterizer.cc \
  CylinderGeom_MVTX.cc

libmvtx_la_LIBADD = \
  libmvtx_io.la \
  -lfun4all \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lg4detectors

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
