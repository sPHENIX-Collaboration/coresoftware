##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_reco_io.la \
  libtrack_reco.la

AM_CPPFLAGS = \
  -DRAVE -DRaveDllExport= \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  ActsAlignmentStates.h \
  ALICEKF.h \
  AssocInfoContainer.h \
  AssocInfoContainerv1.h \
  GPUTPCBaseTrackParam.h \
  GPUTPCTrackLinearisation.h \
  GPUTPCTrackParam.h \
  MakeActsGeometry.h \
  nanoflann.hpp \
  PHActsGSF.h \
  PHActsKDTreeSeeding.h \
  PHActsSiliconSeeding.h \
  PHActsVertexPropagator.h \
  PHActsTrkFitter.h \
  PHActsTrackProjection.h \
  PHCASeeding.h \
  PHGenFitTrackProjection.h \
  PHGenFitTrkFitter.h \
  PHGhostRejection.h \
  PHInitVertexing.h \
  PHMicromegasTpcTrackMatching.h \
  PHSiliconTpcTrackMatching.h \
  PHSimpleVertexFinder.h \
  PHTrackCleaner.h \
  PHTrackSelector.h \
  PHRaveVertexing.h \
  PHSiliconSeedMerger.h \
  PHSiliconTruthTrackSeeding.h \
  PHSimpleKFProp.h \
  PHTpcClusterMover.h \
  PHTpcDeltaZCorrection.h \
  PHTrackClusterAssociator.h \
  PHTrackPropagating.h \
  PHTrackSeeding.h \
  PHTrackFitting.h \
  PHTrackSetMerging.h \
  PHTrackSetCopyMerging.h \
  PHTruthClustering.h \
  PHTruthTrackSeeding.h \
  PHTruthSiliconAssociation.h \
  PHTruthVertexing.h \
  SvtxTrackStateRemoval.h \
  TrackingIterationCounter.h

ROOTDICTS = \
  AssocInfoContainer_Dict.cc \
  AssocInfoContainerv1_Dict.cc

  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  AssocInfoContainer_Dict_rdict.pcm \
  AssocInfoContainerv1_Dict_rdict.pcm

libtrack_reco_io_la_SOURCES = \
  $(ROOTDICTS) \
  AssocInfoContainerv1.cc \
  AssocInfoContainer.cc \
  GPUTPCTrackParam.cxx

ACTS_SOURCES = \
  ActsAlignmentStates.cc \
  MakeActsGeometry.cc \
  PHActsGSF.cc \
  PHActsKDTreeSeeding.cc \
  PHActsSiliconSeeding.cc \
  PHActsTrkFitter.cc \
  PHActsVertexPropagator.cc \
  PHActsTrackProjection.cc

$OFFLINE_MAIN/share:
dist_data_DATA = \
  tgeo-sphenix.json \
  tgeo-sphenix-mms.json

AM_CPPFLAGS += -I$(OFFLINE_MAIN)/include/ActsFatras

ACTS_LIBS = \
  -lActsCore \
  -lActsPluginTGeo \
  -lActsExamplesCommon \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework

libtrack_reco_la_SOURCES = \
  $(ACTS_SOURCES) \
  ALICEKF.cc \
  PH3DVertexing.cc \
  PHCASeeding.cc \
  PHGenFitTrackProjection.cc \
  PHGenFitTrkFitter.cc \
  PHGhostRejection.cc \
  PHInitVertexing.cc \
  PHMicromegasTpcTrackMatching.cc \
  PHRaveVertexing.cc \
  PHSiliconSeedMerger.cc \
  PHSiliconTpcTrackMatching.cc \
  PHSiliconTruthTrackSeeding.cc \
  PHSimpleKFProp.cc \
  PHSimpleVertexFinder.cc \
  PHTpcClusterMover.cc \
  PHTpcDeltaZCorrection.cc \
  PHTrackCleaner.cc \
  PHTrackClusterAssociator.cc \
  PHTrackSeeding.cc \
  PHTrackSelector.cc \
  PHTrackSetMerging.cc \
  PHTrackPropagating.cc \
  PHTrackFitting.cc \
  PHTruthClustering.cc \
  PHTruthTrackSeeding.cc \
  PHTruthVertexing.cc \
  PHTruthSiliconAssociation.cc \
  SvtxTrackStateRemoval.cc \
  TrackingIterationCounter.cc

libtrack_reco_io_la_LIBADD = \
  -lphool

libtrack_reco_la_LIBADD = \
  libtrack_reco_io.la \
  $(ACTS_LIBS) \
  -lSubsysReco \
  -lg4eval \
  -lg4testbench \
  -lg4detectors \
  -lgsl \
  -lgslcblas \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
  -lg4bbc_io \
  -lg4tpc \
  -lg4intt \
  -lg4mvtx \
  -lmicromegas_io \
  -lmvtx_io \
  -lintt_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lcalo_io \
  -lphparameter


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_track_reco \
  testexternals_track_reco_io


testexternals_track_reco_SOURCES = testexternals.cc
testexternals_track_reco_LDADD = libtrack_reco.la

testexternals_track_reco_io_SOURCES = testexternals.cc
testexternals_track_reco_io_LDADD = libtrack_reco_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
