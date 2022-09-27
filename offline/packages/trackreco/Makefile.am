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
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  ../PHTpcTracker/PHTpcTrackerUtil.h \
  ../PHTpcTracker/externals/kdfinder.hpp \
  ActsEvaluator.h \
  ALICEKF.h \
  AssocInfoContainer.h \
  AssocInfoContainerv1.h \
  CellularAutomaton.h \
  CellularAutomaton_v1.h \
  GPUTPCBaseTrackParam.h \
  GPUTPCTrackLinearisation.h \
  GPUTPCTrackParam.h \
  HelixHoughBin.h \
  HelixHoughBin_v1.h \
  HelixHoughFuncs.h \
  HelixHoughFuncs_v1.h \
  HelixHoughSpace.h \
  HelixHoughSpace_v1.h \
  HelixKalmanFilter.h \
  MakeActsGeometry.h \
  nanoflann.hpp \
  PHActsGSF.h \
  PHActsSiliconSeeding.h \
  PHActsVertexPropagator.h \
  PHActsTrkFitter.h \
  PHActsTrackProjection.h \
  PHCASeeding.h \
  PHGenFitTrackProjection.h \
  PHGenFitTrkFitter.h \
  PHGhostRejection.h \
  PHInitVertexing.h \
  PHInitZVertexing.h \
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
  VertexFitter.h

ROOTDICTS = \
  AssocInfoContainer_Dict.cc \
  AssocInfoContainerv1_Dict.cc \
  HelixHoughBin_Dict.cc \
  HelixHoughBin_v1_Dict.cc \
  HelixHoughFuncs_Dict.cc \
  HelixHoughFuncs_v1_Dict.cc \
  HelixHoughSpace_Dict.cc \
  HelixHoughSpace_v1_Dict.cc

  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  AssocInfoContainer_Dict_rdict.pcm \
  AssocInfoContainerv1_Dict_rdict.pcm \
  HelixHoughBin_Dict_rdict.pcm \
  HelixHoughBin_v1_Dict_rdict.pcm \
  HelixHoughFuncs_Dict_rdict.pcm \
  HelixHoughFuncs_v1_Dict_rdict.pcm \
  HelixHoughSpace_Dict_rdict.pcm \
  HelixHoughSpace_v1_Dict_rdict.pcm

libtrack_reco_io_la_SOURCES = \
  $(ROOTDICTS) \
  AssocInfoContainerv1.cc \
  AssocInfoContainer.cc \
  CellularAutomaton_v1.cc \
  GPUTPCTrackParam.cxx \
  HelixHoughBin.cc \
  HelixHoughBin_v1.cc \
  HelixHoughFuncs_v1.cc \
  HelixHoughSpace_v1.cc \
  HelixKalmanFilter.cc \
  VertexFitter.cc

if MAKE_ACTS
ACTS_SOURCES = \
  ActsEvaluator.cc \
  MakeActsGeometry.cc \
  PHActsGSF.cc \
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
  -lActsFatras \
  -lActsPluginTGeo \
  -lActsExamplesCommon \
  -lActsPluginJson \
  -lActsExamplesTrackFitting \
  -lActsExamplesTrackFinding \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework \
  -lActsExamplesIoRoot \
  -lActsExamplesMagneticField
endif

libtrack_reco_la_SOURCES = \
  $(ACTS_SOURCES) \
  ../PHTpcTracker/PHTpcTrackerUtil.cc \
  ALICEKF.cc \
  PH3DVertexing.cc \
  PHCASeeding.cc \
  PHGenFitTrackProjection.cc \
  PHGenFitTrkFitter.cc \
  PHGhostRejection.cc \
  PHInitVertexing.cc \
  PHInitZVertexing.cc \
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
  PHTrackSeeding.cc \
  PHTrackSelector.cc \
  PHTrackSetMerging.cc \
  PHTrackPropagating.cc \
  PHTrackFitting.cc \
  PHTruthClustering.cc \
  PHTruthTrackSeeding.cc \
  PHTruthVertexing.cc \
  PHTruthSiliconAssociation.cc \
  SvtxTrackStateRemoval.cc

libtrack_reco_io_la_LIBADD = \
  -lphool \
  -lHelixHough 

libtrack_reco_la_LIBADD = \
  libtrack_reco_io.la \
  $(ACTS_LIBS) \
  -lSubsysReco \
  -lg4eval \
  -lg4testbench \
  -lg4detectors \
  -lgsl \
  -lgslcblas \
  -lFitNewton \
  -lSeamstress \
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
  -ltrackbase_historic_io \
  -lcalo_io \
  -lphparameter \
  -llog4cpp


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
