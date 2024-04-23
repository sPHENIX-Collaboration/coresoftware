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
  TrainingHitsContainer.h \
  TrainingHits.h \
  TpcRawDataTree.h \
  TpcClusterCleaner.h \
  TpcClusterizer.h \
  TpcClusterMover.h \
  TpcClusterZCrossingCorrection.h \
  TpcCombinedRawDataUnpacker.h \
  TpcDistortionCorrection.h \
  TpcDistortionCorrectionContainer.h \
  TpcLoadDistortionCorrection.h \
  TpcMap.h \
  TpcRawWriter.h \
  TpcSimpleClusterizer.h

ROOTDICTS = \
  LaserClusterizer_Dict.cc \
  TrainingHitsContainer_Dict.cc \
  TrainingHits_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  LaserClusterizer_Dict_rdict.pcm \
  TrainingHitsContainer_Dict_rdict.pcm \
  TrainingHits_Dict_rdict.pcm

mydatadir = $(datadir)/$(PACKAGE)
dist_mydata_DATA = \
  net_model.pt

# sources for tpc library
libtpc_la_SOURCES = \
  LaserClusterizer.cc \
  TpcRawDataTree.cc \
  TpcClusterCleaner.cc \
  TpcClusterizer.cc \
  TpcCombinedRawDataUnpacker.cc \
  TpcLoadDistortionCorrection.cc \
  TpcMap.cc \
  TpcRawWriter.cc \
  TpcSimpleClusterizer.cc

libtpc_la_LIBADD = \
  libtpc_io.la \
  -lffarawobjects \
  -lc10 \
  -ltorch_cpu \
  -ltorch \
  -lmicromegas_io \
  -lfun4all \
  -lSpectrum \
  -lpthread

# sources for io library
libtpc_io_la_SOURCES = \
  $(ROOTDICTS) \
  LaserClusterizer.cc \
  TrainingHitsContainer.cc \
  TrainingHits.cc \
  TpcClusterMover.cc \
  TpcDistortionCorrection.cc \
  TpcClusterZCrossingCorrection.cc

libtpc_io_la_LIBADD = \
  -ltrack_io \
  -lg4detectors_io \
  -ltrackbase_historic_io

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
