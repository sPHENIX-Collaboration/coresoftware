##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmvtx_io.la \
  libmvtx.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  CylinderGeom_Mvtx.h \
  CylinderGeom_MvtxHelper.h \
  MvtxCombinedRawDataDecoder.h \
  MvtxClusterizer.h \
  MvtxHitPruner.h \
  MvtxNoiseMap.h \
  MvtxPixelDefs.h \
  MvtxPixelMask.h \
  SegmentationAlpide.h

ROOTDICTS = \
  CylinderGeom_Mvtx_Dict.cc \
  MvtxNoiseMap_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CylinderGeom_Mvtx_Dict_rdict.pcm \
  MvtxNoiseMap_Dict_rdict.pcm

# sources for mvtx library
libmvtx_la_SOURCES = \
  CylinderGeom_MvtxHelper.cc \
  MvtxCombinedRawDataDecoder.cc \
  MvtxClusterizer.cc \
  MvtxHitPruner.cc \
  MvtxPixelDefs.cc \
  MvtxPixelMask.cc

libmvtx_la_LIBADD = \
  libmvtx_io.la \
  -lActsCore \
  -lCLHEP \
  -lcdbobjects \
  -lffamodules \
  -lffarawobjects \
  -lfun4all \
  -lfun4allraw \
  -lmvtx_decoder \
  -lphg4hit \
  -lSubsysReco \
  -ltrack_io \
  -ltrackbase_historic_io

# sources for io library
libmvtx_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeom_Mvtx.cc \
  MvtxNoiseMap.cc \
  SegmentationAlpide.cc

libmvtx_io_la_LIBADD = \
  -lg4detectors_io \
  -lphool

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

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
