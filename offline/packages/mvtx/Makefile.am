##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
  libmvtx_decoder.la \
  libmvtx_io.la \
  libmvtx.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  MvtxClusterizer.h \
  MvtxClusterQA.h \
  MvtxHitPruner.h \
  MvtxCombinedRawDataDecoder.h \
  CylinderGeom_Mvtx.h \
  SegmentationAlpide.h

decoderincludedir = $(includedir)/mvtx_decoder
decoderinclude_HEADERS = \
  mvtx_decoder/DecodingStat.h \
  mvtx_decoder/GBTLink.h \
  mvtx_decoder/GBTWord.h \
  mvtx_decoder/InteractionRecord.h \
  mvtx_decoder/mvtx_utils.h \
  mvtx_decoder/PayLoadCont.h \
  mvtx_decoder/PayLoadSG.h \
  mvtx_decoder/PixelData.h

ROOTDICTS = \
  CylinderGeom_Mvtx_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CylinderGeom_Mvtx_Dict_rdict.pcm

# source for mvtx decoder library
libmvtx_decoder_la_SOURCES =\
  mvtx_decoder/InteractionRecord.cc \
  mvtx_decoder/PayLoadCont.cc \
  mvtx_decoder/PixelData.cc \
  mvtx_decoder/GBTWord.cc \
  mvtx_decoder/DecodingStat.cc \
  mvtx_decoder/GBTLink.cc


# sources for mvtx library
libmvtx_la_SOURCES = \
  MvtxClusterizer.cc \
  MvtxClusterQA.cc \
  MvtxHitPruner.cc \
  MvtxCombinedRawDataDecoder.cc

libmvtx_la_LIBADD = \
  libmvtx_decoder.la \
  libmvtx_io.la \
  -lCLHEP \
  -lffarawobjects \
  -lphg4hit \
  -lSubsysReco \
  -lqautils \
  -lcdbobjects

# sources for io library
libmvtx_io_la_SOURCES = \
  $(ROOTDICTS) \
  CylinderGeom_Mvtx.cc \
  SegmentationAlpide.cc

libmvtx_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lcdbobjects

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;



################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_mvtx_decoder \
  testexternals_mvtx_io \
  testexternals_mvtx

testexternals_mvtx_decoder_SOURCES = testexternals.cc
testexternals_mvtx_decoder_LDADD = libmvtx_decoder.la

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
