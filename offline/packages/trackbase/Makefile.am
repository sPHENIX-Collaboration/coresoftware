##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  ActsGeometry.h \
  ActsGsfTrackFittingAlgorithm.h \
  ActsSourceLink.h \
  ActsSurfaceMaps.h \
  ActsTrackFittingAlgorithm.h \
  ActsTrackingGeometry.h \
  AlignmentTransformation.h \
  CMFlashCluster.h \
  CMFlashClusterContainer.h \
  CMFlashClusterContainerv1.h \
  CMFlashClusterv1.h \
  CMFlashDifference.h \
  CMFlashDifferenceContainer.h \
  CMFlashDifferenceContainerv1.h \
  CMFlashDifferencev1.h \
  Calibrator.h \
  ClusterErrorPara.h \
  EmbRecoMatch.h \
  EmbRecoMatchv1.h \
  EmbRecoMatchContainer.h \
  EmbRecoMatchContainerv1.h \
  InttDefs.h \
  MvtxDefs.h \
  RawHit.h \
  RawHitSet.h \
  RawHitSetContainer.h \
  RawHitSetContainerv1.h \
  RawHitSetv1.h \
  RawHitTpc.h \
  RawHitv1.h \
  ResidualOutlierFinder.h \
  SpacePoint.h \
  TpcDefs.h \
  TpcSeedTrackMap.h \
  TpcSeedTrackMapv1.h \
  TrackFitUtils.h \
  TrkrCluster.h \
  TrkrClusterContainer.h \
  TrkrClusterContainerv1.h \
  TrkrClusterContainerv2.h \
  TrkrClusterContainerv3.h \
  TrkrClusterContainerv4.h \
  TrkrClusterCrossingAssoc.h \
  TrkrClusterCrossingAssocv1.h \
  TrkrClusterHitAssoc.h \
  TrkrClusterHitAssocv1.h \
  TrkrClusterHitAssocv2.h \
  TrkrClusterHitAssocv3.h \
  TrkrClusterIterationMap.h \
  TrkrClusterIterationMapv1.h \
  TrkrClusterv1.h \
  TrkrClusterv2.h \
  TrkrClusterv3.h \
  TrkrClusterv4.h \
  TrkrDefs.h \
  TrkrHit.h \
  TrkrHitSet.h \
  TrkrHitSetContainer.h \
  TrkrHitSetContainerv1.h \
  TrkrHitSetv1.h \
  TrkrHitTruthAssoc.h \
  TrkrHitTruthAssocv1.h \
  TrkrHitTruthClusters.h \
  TrkrHitTruthClustersv1.h \
  TrkrHitv1.h \
  TrkrHitv2.h \
  TrkrTruthClusters.h \
  TrkrTruthClustersv1.h \
  TrkrTruthTrack.h \
  TrkrTruthTrackContainer.h \
  TrkrTruthTrackContainerv1.h \
  TrkrTruthTrackv1.h \
  alignmentTransformationContainer.h \
  sPHENIXActsDetectorElement.h

ROOTDICTS = \
  CMFlashClusterContainer_Dict.cc \
  CMFlashClusterContainerv1_Dict.cc \
  CMFlashCluster_Dict.cc \
  CMFlashClusterv1_Dict.cc \
  CMFlashDifferenceContainer_Dict.cc \
  CMFlashDifferenceContainerv1_Dict.cc \
  CMFlashDifference_Dict.cc \
  CMFlashDifferencev1_Dict.cc \
  EmbRecoMatchContainer_Dict.cc \
  EmbRecoMatchContainerv1_Dict.cc \
  EmbRecoMatch_Dict.cc \
  EmbRecoMatchv1_Dict.cc \
  RawHitSetContainer_Dict.cc \
  RawHitSetContainerv1_Dict.cc \
  RawHitSet_Dict.cc \
  RawHitSetv1_Dict.cc \
  RawHitTpc_Dict.cc \
  RawHit_Dict.cc \
  RawHitv1_Dict.cc \
  TpcSeedTrackMap_Dict.cc \
  TpcSeedTrackMapv1_Dict.cc \
  TrkrClusterContainer_Dict.cc \
  TrkrClusterContainerv1_Dict.cc \
  TrkrClusterContainerv2_Dict.cc \
  TrkrClusterContainerv3_Dict.cc \
  TrkrClusterContainerv4_Dict.cc \
  TrkrClusterCrossingAssoc_Dict.cc \
  TrkrClusterCrossingAssocv1_Dict.cc \
  TrkrClusterHitAssoc_Dict.cc \
  TrkrClusterHitAssocv1_Dict.cc \
  TrkrClusterHitAssocv2_Dict.cc \
  TrkrClusterHitAssocv3_Dict.cc \
  TrkrClusterIterationMap_Dict.cc \
  TrkrClusterIterationMapv1_Dict.cc \
  TrkrCluster_Dict.cc \
  TrkrClusterv1_Dict.cc \
  TrkrClusterv2_Dict.cc \
  TrkrClusterv3_Dict.cc \
  TrkrClusterv4_Dict.cc \
  TrkrHitSetContainer_Dict.cc \
  TrkrHitSetContainerv1_Dict.cc \
  TrkrHitSet_Dict.cc \
  TrkrHitSetv1_Dict.cc \
  TrkrHitTruthAssoc_Dict.cc \
  TrkrHitTruthAssocv1_Dict.cc \
  TrkrHitTruthClusters_Dict.cc \
  TrkrHitTruthClustersv1_Dict.cc \
  TrkrHit_Dict.cc \
  TrkrHitv1_Dict.cc \
  TrkrHitv2_Dict.cc \
  TrkrTruthClusters_Dict.cc \
  TrkrTruthClustersv1_Dict.cc \
  TrkrTruthTrackContainer_Dict.cc \
  TrkrTruthTrackContainerv1_Dict.cc \
  TrkrTruthTrack_Dict.cc \
  TrkrTruthTrackv1_Dict.cc \
  alignmentTransformationContainer_Dict.cc \
  sPHENIXActsDetectorElement.cc


pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CMFlashClusterContainer_Dict_rdict.pcm \
  CMFlashClusterContainerv1_Dict_rdict.pcm \
  CMFlashCluster_Dict_rdict.pcm \
  CMFlashClusterv1_Dict_rdict.pcm \
  CMFlashDifferenceContainer_Dict_rdict.pcm \
  CMFlashDifferenceContainerv1_Dict_rdict.pcm \
  CMFlashDifference_Dict_rdict.pcm \
  CMFlashDifferencev1_Dict_rdict.pcm \
  EmbRecoMatchContainer_Dict_rdict.pcm \
  EmbRecoMatchContainerv1_Dict_rdict.pcm \
  EmbRecoMatch_Dict_rdict.pcm \
  EmbRecoMatchv1_Dict_rdict.pcm \
  RawHitSetContainer_Dict_rdict.pcm \
  RawHitSetContainerv1_Dict_rdict.pcm \
  RawHitSet_Dict_rdict.pcm \
  RawHitSetv1_Dict_rdict.pcm \
  RawHitTpc_Dict_rdict.pcm \
  RawHit_Dict_rdict.pcm \
  RawHitv1_Dict_rdict.pcm \
  TpcSeedTrackMap_Dict_rdict.pcm \
  TpcSeedTrackMapv1_Dict_rdict.pcm \
  TrkrClusterContainer_Dict_rdict.pcm \
  TrkrClusterContainerv1_Dict_rdict.pcm \
  TrkrClusterContainerv2_Dict_rdict.pcm \
  TrkrClusterContainerv3_Dict_rdict.pcm \
  TrkrClusterContainerv4_Dict_rdict.pcm \
  TrkrClusterCrossingAssoc_Dict_rdict.pcm \
  TrkrClusterCrossingAssocv1_Dict_rdict.pcm \
  TrkrClusterHitAssoc_Dict_rdict.pcm \
  TrkrClusterHitAssocv1_Dict_rdict.pcm \
  TrkrClusterHitAssocv2_Dict_rdict.pcm \
  TrkrClusterHitAssocv3_Dict_rdict.pcm \
  TrkrClusterIterationMap_Dict_rdict.pcm \
  TrkrClusterIterationMapv1_Dict_rdict.pcm \
  TrkrCluster_Dict_rdict.pcm \
  TrkrClusterv1_Dict_rdict.pcm \
  TrkrClusterv2_Dict_rdict.pcm \
  TrkrClusterv3_Dict_rdict.pcm \
  TrkrClusterv4_Dict_rdict.pcm \
  TrkrHitSetContainer_Dict_rdict.pcm \
  TrkrHitSetContainerv1_Dict_rdict.pcm \
  TrkrHitSet_Dict_rdict.pcm \
  TrkrHitSetv1_Dict_rdict.pcm \
  TrkrHitTruthAssoc_Dict_rdict.pcm \
  TrkrHitTruthAssocv1_Dict_rdict.pcm \
  TrkrHitTruthClusters_Dict_rdict.pcm \
  TrkrHitTruthClustersv1_Dict_rdict.pcm \
  TrkrHit_Dict_rdict.pcm \
  TrkrHitv1_Dict_rdict.pcm \
  TrkrHitv2_Dict_rdict.pcm \
  TrkrTruthClusters_Dict_rdict.pcm \
  TrkrTruthClustersv1_Dict_rdict.pcm \
  TrkrTruthTrackContainer_Dict_rdict.pcm \
  TrkrTruthTrackContainerv1_Dict_rdict.pcm \
  TrkrTruthTrack_Dict_rdict.pcm \
  TrkrTruthTrackv1_Dict_rdict.pcm \
  alignmentTransformationContainer_Dict_rdict.pcm

# sources for io library
libtrack_io_la_SOURCES = \
  $(ROOTDICTS) \
  ActsGeometry.cc \
  ActsSurfaceMaps.cc \
  AlignmentTransformation.cc \
  CMFlashClusterContainerv1.cc \
  CMFlashClusterv1.cc \
  CMFlashDifferenceContainerv1.cc \
  CMFlashDifferencev1.cc \
  ClusterErrorPara.cc \
  EmbRecoMatchContainer.cc \
  EmbRecoMatchContainerv1.cc \
  EmbRecoMatchv1.cc \
  InttDefs.cc \
  MvtxDefs.cc \
  RawHitSet.cc \
  RawHitSetContainer.cc \
  RawHitSetContainerv1.cc \
  RawHitSetv1.cc \
  RawHitTpc.cc \
  RawHitv1.cc \
  TpcDefs.cc \
  TpcSeedTrackMap.cc \
  TpcSeedTrackMapv1.cc \
  TrackFitUtils.cc \
  TrackFittingAlgorithmFunctionsGsf.cc \
  TrackFittingAlgorithmFunctionsKalman.cc \
  TrkrClusterContainer.cc \
  TrkrClusterContainerv1.cc \
  TrkrClusterContainerv2.cc \
  TrkrClusterContainerv3.cc \
  TrkrClusterContainerv4.cc \
  TrkrClusterCrossingAssoc.cc \
  TrkrClusterCrossingAssocv1.cc \
  TrkrClusterHitAssoc.cc \
  TrkrClusterHitAssocv1.cc \
  TrkrClusterHitAssocv2.cc \
  TrkrClusterHitAssocv3.cc \
  TrkrClusterIterationMap.cc \
  TrkrClusterIterationMapv1.cc \
  TrkrClusterv1.cc \
  TrkrClusterv2.cc \
  TrkrClusterv3.cc \
  TrkrClusterv4.cc \
  TrkrDefs.cc \
  TrkrHitSet.cc \
  TrkrHitSetContainer.cc \
  TrkrHitSetContainerv1.cc \
  TrkrHitSetv1.cc \
  TrkrHitTruthAssocv1.cc \
  TrkrHitTruthClustersv1.cc \
  TrkrHitv1.cc \
  TrkrHitv2.cc \
  TrkrTruthClusters.cc \
  TrkrTruthClustersv1.cc \
  TrkrTruthTrack.cc \
  TrkrTruthTrackContainer.cc \
  TrkrTruthTrackContainerv1.cc \
  TrkrTruthTrackv1.cc \
  alignmentTransformationContainer.cc

libtrack_io_la_LIBADD = \
  -lphool \
  -lActsCore \
  -lActsPluginTGeo \
  -lffamodules


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_track_io

#testexternals_track_util_SOURCES = testexternals.cc
#testexternals_track_util_LDADD = libtrack_util.la

testexternals_track_io_SOURCES = testexternals.cc
testexternals_track_io_LDADD = libtrack_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
