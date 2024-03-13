AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  Fun4AllEventOutStream.h \
  Fun4AllEventOutputManager.h \
  Fun4AllFileOutStream.h \
  Fun4AllPrdfInputManager.h \
  Fun4AllPrdfInputPoolManager.h \
  Fun4AllPrdfInputTriggerManager.h \
  Fun4AllPrdfOutputManager.h \
  Fun4AllRolloverFileOutStream.h \
  Fun4AllStreamingInputManager.h \
  InputManagerType.h \
  SingleGl1PoolInput.h \
  SingleGl1TriggerInput.h \
  SingleHcalTriggerInput.h \
  SingleMbdTriggerInput.h \
  SingleInttPoolInput.h \
  SingleMicromegasPoolInput.h \
  SingleMvtxPoolInput.h \
  SinglePrdfInput.h \
  SingleStreamingInput.h \
  SingleTpcPoolInput.h \
  SingleTriggerInput.h \
  SingleZdcInput.h

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

lib_LTLIBRARIES = \
  libmvtx_decoder.la \
  libfun4allraw.la

# source for mvtx decoder library
libmvtx_decoder_la_SOURCES =\
  mvtx_decoder/InteractionRecord.cc \
  mvtx_decoder/PayLoadCont.cc \
  mvtx_decoder/PixelData.cc \
  mvtx_decoder/GBTWord.cc \
  mvtx_decoder/DecodingStat.cc \
  mvtx_decoder/GBTLink.cc

libfun4allraw_la_SOURCES = \
  Fun4AllEventOutStream.cc \
  Fun4AllEventOutputManager.cc \
  Fun4AllFileOutStream.cc \
  Fun4AllPrdfInputManager.cc \
  Fun4AllPrdfInputPoolManager.cc \
  Fun4AllPrdfInputTriggerManager.cc \
  Fun4AllPrdfOutputManager.cc \
  Fun4AllRolloverFileOutStream.cc \
  Fun4AllStreamingInputManager.cc \
  intt_pool.cc \
  mvtx_pool.cc \
  tpc_pool.cc \
  SingleGl1PoolInput.cc \
  SingleGl1TriggerInput.cc \
  SingleHcalTriggerInput.cc \
  SingleInttPoolInput.cc \
  SingleMbdTriggerInput.cc \
  SingleMicromegasPoolInput.cc \
  SingleMvtxPoolInput.cc \
  SinglePrdfInput.cc \
  SingleStreamingInput.cc \
  SingleTpcPoolInput.cc \
  SingleTriggerInput.cc \
  SingleZdcInput.cc

libfun4allraw_la_LIBADD = \
  libmvtx_decoder.la \
  -lffarawobjects \
  -lfun4all \
  -lEvent \
  -lphoolraw

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_mvtx_decoder \
  testexternals

testexternals_mvtx_decoder_SOURCES = testexternals.cc
testexternals_mvtx_decoder_LDADD = libmvtx_decoder.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libfun4allraw.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
