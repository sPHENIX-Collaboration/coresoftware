##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

if USE_ONLINE

lib_LTLIBRARIES = \
  libtpc.la

pkginclude_HEADERS = \
  TpcMap.h

libtpc_la_SOURCES = \
  TpcMap.cc

noinst_PROGRAMS = \
  testexternals_tpc

else

lib_LTLIBRARIES = \
  libtpc_io.la \
  libtpc.la

pkginclude_HEADERS = \
  LaserClusterizer.h \
  LaserEventInfo.h \
  LaserEventInfov1.h \
  LaserEventIdentifier.h \
  LaserEventRejecter.h \
  TrainingHitsContainer.h \
  TrainingHits.h \
  Tpc3DClusterizer.h \
  TpcRawDataTree.h \
  TpcClusterCleaner.h \
  TpcClusterizer.h \
  TpcClusterMover.h \
  TpcClusterZCrossingCorrection.h \
  TpcCombinedRawDataUnpacker.h \
  TpcDistortionCorrection.h \
  TpcDistortionCorrectionContainer.h \
  TpcGlobalPositionWrapper.h \
  TpcLoadDistortionCorrection.h \
  TpcMap.h \
  TpcRawWriter.h \
  TpcSimpleClusterizer.h

ROOTDICTS = \
  LaserEventInfo_Dict.cc \
  LaserEventInfov1_Dict.cc \
  TrainingHitsContainer_Dict.cc \
  TrainingHits_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  LaserEventInfo_Dict_rdict.pcm \
  LaserEventInfov1_Dict_rdict.pcm \
  TrainingHitsContainer_Dict_rdict.pcm \
  TrainingHits_Dict_rdict.pcm

mydatadir = $(datadir)/$(PACKAGE)
dist_mydata_DATA = \
  net_model.pt

# sources for tpc library
libtpc_la_SOURCES = \
  LaserClusterizer.cc \
  LaserEventInfov1.cc \
  LaserEventIdentifier.cc \
  LaserEventRejecter.cc \
  TpcRawDataTree.cc \
  Tpc3DClusterizer.cc \
  TpcClusterCleaner.cc \
  TpcClusterizer.cc \
  TpcCombinedRawDataUnpacker.cc \
  TpcGlobalPositionWrapper.cc \
  TpcLoadDistortionCorrection.cc \
  TpcMap.cc \
  TpcRawWriter.cc \
  TpcSimpleClusterizer.cc \
  TpcClusterMover.cc \
  TpcClusterZCrossingCorrection.cc \
  TpcDistortionCorrection.cc

libtpc_la_LIBADD = \
  libtpc_io.la \
  -lffarawobjects \
  -lc10 \
  -ltorch_cpu \
  -ltorch \
  -ltrack \
  -lmicromegas_io \
  -lfun4all \
  -lphool \
  -lSpectrum \
  -lpthread

# sources for io library
libtpc_io_la_SOURCES = \
  $(ROOTDICTS) \
  LaserEventInfov1.cc \
  TrainingHitsContainer.cc \
  TrainingHits.cc

libtpc_io_la_LIBADD = \
  -lg4detectors_io

noinst_PROGRAMS = \
  testexternals_tpc_io \
  testexternals_tpc

endif

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

testexternals_tpc_io_SOURCES = testexternals.cc
testexternals_tpc_io_LDADD = libtpc_io.la

testexternals_tpc_SOURCES = testexternals.cc
testexternals_tpc_LDADD = libtpc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* *.pcm $(BUILT_SOURCES)
