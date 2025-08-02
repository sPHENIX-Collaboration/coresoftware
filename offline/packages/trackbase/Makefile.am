##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign


if USE_ONLINE
# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_io.la

pkginclude_HEADERS = \
  TrkrDefs.h

libtrack_io_la_SOURCES = \
  TrkrDefs.cc

noinst_PROGRAMS = \
  testexternals_track_io

else

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_io.la \
  libtrack.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib


pkginclude_HEADERS = \
  ActsAborter.h \
  ActsGeometry.h \
  ActsGsfTrackFittingAlgorithm.h \
  ActsSourceLink.h \
  ActsSurfaceMaps.h \
  ActsTrackFittingAlgorithm.h \
  ActsTrackingGeometry.h \
  alignmentTransformationContainer.h \
  AlignmentTransformation.h \
  CommonOptions.h \
  CMFlashCluster.h \
  CMFlashClusterContainer.h \
  CMFlashClusterContainerv1.h \
  CMFlashClusterv1.h \
  CMFlashClusterv2.h \
  CMFlashClusterv3.h \
  CMFlashDifference.h \
  CMFlashDifferenceContainer.h \
  CMFlashDifferenceContainerv1.h \
  CMFlashDifferencev1.h \
  Calibrator.h \
  ClusHitsVerbose.h \
  ClusHitsVerbosev1.h \
  ClusterErrorPara.h \
  IBaseDetector.h \
  InttDefs.h \
  InttEventInfo.h \
  InttEventInfov1.h \
  LaserCluster.h \
  LaserClusterContainer.h \
  LaserClusterContainerv1.h \
  LaserClusterv1.h \
  LaserClusterv2.h \
  MagneticFieldOptions.h \
  MaterialWiper.h \
  MvtxDefs.h \
  MvtxEventInfo.h \
  MvtxEventInfov1.h \
  MvtxEventInfov2.h \
  MvtxEventInfov3.h \
  RawHit.h \
  RawHitSet.h \
  RawHitSetContainer.h \
  RawHitSetContainerv1.h \
  RawHitSetv1.h \
  RawHitTpc.h \
  RawHitv1.h \
  ResidualOutlierFinder.h \
  SpacePoint.h \
  sPHENIXActsDetectorElement.h \
  TGeoDetectorWithOptions.h \
  TpcDefs.h \
  TpcSeedTrackMap.h \
  TpcSeedTrackMapv1.h \
  TpcTpotEventInfo.h \
  TpcTpotEventInfov1.h \
  TrackVertexCrossingAssoc.h \
  TrackVertexCrossingAssoc_v1.h \
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
  TrkrClusterv5.h \
  TrkrDefs.h \
  TrkrHit.h \
  TrkrHitSet.h \
  TrkrHitSetContMvtxHelper.h \
  TrkrHitSetContMvtxHelperv1.h \
  TrkrHitSetContainer.h \
  TrkrHitSetContainerv1.h \
  TrkrHitSetContainerv2.h \
  TrkrHitSetv1.h \
  TrkrHitSetTpc.h \
  TrkrHitSetTpcv1.h \
  TrkrHitTruthAssoc.h \
  TrkrHitTruthAssocv1.h \
  TrkrHitv1.h \
  TrkrHitv2.h

ROOTDICTS = \
  CMFlashClusterContainer_Dict.cc \
  CMFlashClusterContainerv1_Dict.cc \
  CMFlashCluster_Dict.cc \
  CMFlashClusterv1_Dict.cc \
  CMFlashClusterv2_Dict.cc \
  CMFlashClusterv3_Dict.cc \
  CMFlashDifferenceContainer_Dict.cc \
  CMFlashDifferenceContainerv1_Dict.cc \
  CMFlashDifference_Dict.cc \
  CMFlashDifferencev1_Dict.cc \
  ClusHitsVerbose_Dict.cc \
  ClusHitsVerbosev1_Dict.cc \
  InttEventInfo_Dict.cc \
  InttEventInfov1_Dict.cc \
  LaserClusterContainer_Dict.cc \
  LaserClusterContainerv1_Dict.cc \
  LaserCluster_Dict.cc \
  LaserClusterv1_Dict.cc \
  LaserClusterv2_Dict.cc \
  MvtxEventInfo_Dict.cc \
  MvtxEventInfov1_Dict.cc \
  MvtxEventInfov2_Dict.cc \
  MvtxEventInfov3_Dict.cc \
  RawHitSetContainer_Dict.cc \
  RawHitSetContainerv1_Dict.cc \
  RawHitSet_Dict.cc \
  RawHitSetv1_Dict.cc \
  RawHitTpc_Dict.cc \
  RawHit_Dict.cc \
  RawHitv1_Dict.cc \
  TpcSeedTrackMap_Dict.cc \
  TpcSeedTrackMapv1_Dict.cc \
  TpcTpotEventInfo_Dict.cc \
  TpcTpotEventInfov1_Dict.cc \
  TrackVertexCrossingAssoc_Dict.cc \
  TrackVertexCrossingAssoc_v1_Dict.cc \
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
  TrkrClusterv5_Dict.cc \
  TrkrHitSetContMvtxHelper_Dict.cc \
  TrkrHitSetContMvtxHelperv1_Dict.cc \
  TrkrHitSetContainer_Dict.cc \
  TrkrHitSetContainerv1_Dict.cc \
  TrkrHitSetContainerv2_Dict.cc \
  TrkrHitSet_Dict.cc \
  TrkrHitSetv1_Dict.cc \
  TrkrHitSetTpc_Dict.cc \
  TrkrHitSetTpcv1_Dict.cc \
  TrkrHitTruthAssoc_Dict.cc \
  TrkrHitTruthAssocv1_Dict.cc \
  TrkrHit_Dict.cc \
  TrkrHitv1_Dict.cc \
  TrkrHitv2_Dict.cc


pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

# sources for io library
libtrack_la_SOURCES = \
  ActsGeometry.cc \
  ActsSurfaceMaps.cc \
  AlignmentTransformation.cc \
  alignmentTransformationContainer.cc \
  Calibrator.cc \
  ClusterErrorPara.cc \
  CommonOptions.cc \
  MagneticFieldOptions.cc \
  sPHENIXActsDetectorElement.cc \
  TGeoDetectorWithOptions.cc \
  TrackFittingAlgorithmFunctionsGsf.cc \
  TrackFittingAlgorithmFunctionsKalman.cc \
  TrackFitUtils.cc

# sources for io library
libtrack_io_la_SOURCES = \
  $(ROOTDICTS) \
  CMFlashClusterContainerv1.cc \
  CMFlashClusterv1.cc \
  CMFlashClusterv2.cc \
  CMFlashClusterv3.cc \
  CMFlashDifferenceContainerv1.cc \
  CMFlashDifferencev1.cc \
  ClusHitsVerbose.cc \
  ClusHitsVerbosev1.cc \
  InttDefs.cc \
  InttEventInfo.cc \
  InttEventInfov1.cc \
  LaserClusterContainerv1.cc \
  LaserClusterv1.cc \
  LaserClusterv2.cc \
  MvtxDefs.cc \
  MvtxEventInfo.cc \
  MvtxEventInfov1.cc \
  MvtxEventInfov2.cc \
  MvtxEventInfov3.cc \
  RawHitSet.cc \
  RawHitSetContainer.cc \
  RawHitSetContainerv1.cc \
  RawHitSetv1.cc \
  RawHitTpc.cc \
  RawHitv1.cc \
  TpcDefs.cc \
  TpcSeedTrackMap.cc \
  TpcSeedTrackMapv1.cc \
  TpcTpotEventInfov1.cc \
  TrackVertexCrossingAssoc.cc \
  TrackVertexCrossingAssoc_v1.cc \
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
  TrkrClusterv5.cc \
  TrkrDefs.cc \
  TrkrHitSet.cc \
  TrkrHitSetContMvtxHelper.cc \
  TrkrHitSetContMvtxHelperv1.cc \
  TrkrHitSetContainer.cc \
  TrkrHitSetContainerv1.cc \
  TrkrHitSetContainerv2.cc \
  TrkrHitSetv1.cc \
  TrkrHitSetTpc.cc \
  TrkrHitSetTpcv1.cc \
  TrkrHitTruthAssocv1.cc \
  TrkrHitv1.cc \
  TrkrHitv2.cc

libtrack_la_LIBADD = \
  libtrack_io.la \
  -lActsCore \
  -lActsExamplesMagneticField \
  -lActsPluginTGeo \
  -lActsExamplesDetectorTGeo \
  -lffamodules \
  -lboost_program_options

libtrack_io_la_LIBADD = \
  -lphool \
  -lphg4hit

noinst_PROGRAMS = \
  testexternals_track \
  testexternals_track_io

testexternals_track_SOURCES = testexternals.cc
testexternals_track_LDADD = libtrack.la

endif

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc


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
