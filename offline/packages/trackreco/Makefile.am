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
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include/eigen3 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  ActsTransformations.h \
  ActsTrack.h \
  ActsTrackingGeometry.h \
  ActsSurfaceMaps.h \
  ActsEvaluator.h \
  AssocInfoContainer.h \
  CellularAutomaton.h \
  CellularAutomaton_v1.h \
  HelixHoughBin.h \
  HelixHoughBin_v1.h \
  HelixHoughFuncs.h \
  HelixHoughFuncs_v1.h \
  HelixHoughSpace.h \
  HelixHoughSpace_v1.h \
  HelixKalmanFilter.h \
  MakeActsGeometry.h \
  PHGenFitTrackProjection.h \
  PHGenFitTrkFitter.h \
  PHActsSourceLinks.h \
  PHActsTracks.h \
  PHActsTrkProp.h \
  PHActsVertexFinder.h \
  PHActsTrkFitter.h \
  PHActsVertexFitter.h \
  PHGenFitTrkProp.h \
  PHHoughSeeding.h \
  PHRTreeSeeding.h \
  PHCASeeding.h \
  PHInitVertexing.h \
  PHInitZVertexing.h \
  PHSiliconTpcTrackMatching.h \
  PHMicromegasTpcTrackMatching.h \
  PHSiliconTruthTrackSeeding.h \
  PHTrackPropagating.h \
  PHTrackSeeding.h \
  PHTrackFitting.h \
  PHTrackSetMerging.h \
  PHTruthTrackSeeding.h \
  PHTruthSiliconAssociation.h \
  PHTruthVertexing.h \
  VertexFitter.h \
  PHRaveVertexing.h \
  GPUTPCBaseTrackParam.h \
  GPUTPCTrackLinearisation.h \
  GPUTPCTrackParam.h

ROOTDICTS = \
  AssocInfoContainer_Dict.cc \
  HelixHoughBin_Dict.cc \
  HelixHoughBin_v1_Dict.cc \
  HelixHoughFuncs_Dict.cc \
  HelixHoughFuncs_v1_Dict.cc \
  HelixHoughSpace_Dict.cc \
  HelixHoughSpace_v1_Dict.cc

  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  AssocInfoContainer_Dict_rdict.pcm \
  HelixHoughBin_Dict_rdict.pcm \
  HelixHoughBin_v1_Dict_rdict.pcm \
  HelixHoughFuncs_Dict_rdict.pcm \
  HelixHoughFuncs_v1_Dict_rdict.pcm \
  HelixHoughSpace_Dict_rdict.pcm \
  HelixHoughSpace_v1_Dict_rdict.pcm

libtrack_reco_io_la_SOURCES = \
  $(ROOTDICTS) \
  AssocInfoContainer.cc \
  CellularAutomaton_v1.cc \
  HelixHoughBin_v1.cc \
  HelixHoughFuncs_v1.cc \
  HelixHoughSpace_v1.cc \
  HelixKalmanFilter.cc \
  VertexFitter.cc \
  GPUTPCTrackParam.cxx

if MAKE_ACTS
ACTS_SOURCES = \
  ActsTransformations.cc \
  ActsEvaluator.cc \
  MakeActsGeometry.cc \
  PHActsSourceLinks.cc \
  PHActsTracks.cc \
  PHActsTrkProp.cc \
  PHActsTrkFitter.cc \
  PHActsVertexFinder.cc \
  PHActsVertexFitter.cc

$OFFLINE_MAIN/share:
dist_data_DATA = \
  tgeo-sphenix.response \
  tgeo-sphenix-mms.response

AM_CPPFLAGS += -I$(OFFLINE_MAIN)/include/ActsFatras

ACTS_LIBS = \
  -lg4eval \
  -lActsCore \
  -lActsFatras \
  -lActsPluginTGeo \
  -lActsExamplesCommon \
  -lActsPluginJson \
  -lActsExamplesFitting \
  -lActsExamplesTrackFinding \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework \
  -lActsExamplesIoRoot \
  -lActsExamplesMagneticField
endif

libtrack_reco_la_SOURCES = \
  $(ACTS_SOURCES) \
  PHInitVertexing.cc \
  PHInitZVertexing.cc \
  PHTrackSeeding.cc \
  PHTrackSetMerging.cc \
  PHTrackPropagating.cc \
  PHTrackFitting.cc \
  PH3DVertexing.cc \
  PHTruthVertexing.cc \
  PHSiliconTpcTrackMatching.cc \
  PHMicromegasTpcTrackMatching.cc \
  PHTruthTrackSeeding.cc \
  PHSiliconTruthTrackSeeding.cc \
  PHTruthSiliconAssociation.cc \
  PHHoughSeeding.cc \
  PHRTreeSeeding.cc \
  PHCASeeding.cc \
  PHGenFitTrkProp.cc \
  PHGenFitTrkFitter.cc \
  PHGenFitTrackProjection.cc \
  PHRaveVertexing.cc 

libtrack_reco_io_la_LIBADD = \
  -lphool \
  -lHelixHough

libtrack_reco_la_LIBADD = \
  libtrack_reco_io.la \
  $(ACTS_LIBS) \
  -lSubsysReco \
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
