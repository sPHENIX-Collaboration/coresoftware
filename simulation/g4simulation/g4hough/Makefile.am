AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I`root-config --incdir` \
  -I$(OFFLINE_MAIN)/include/eigen3

lib_LTLIBRARIES = \
   libg4hough_io.la \
   libg4hough.la

AM_CXXFLAGS = \
 -Werror -msse2 -DRAVE -DRaveDllExport=

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4hough_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io

libg4hough_la_LDFLAGS = \
  `geant4-config --libs`\
  `root-config --libs`

libg4hough_la_LIBADD = \
  -lg4detectors \
  -lfun4all \
  -lHelixHough \
  -lFitNewton \
  -lSeamstress \
  -lg4bbc_io \
  -lcalo_io -lcalo_util \
  -lgenfit2 \
  -lphfield \
  -lgenfit2exp \
  -lPHGenFit \
  -lSpectrum \
  -lg4tpc \
  libg4hough_io.la

pkginclude_HEADERS = \
  SvtxVertex.h \
  SvtxVertex_v1.h \
  SvtxVertexMap.h \
  SvtxVertexMap_v1.h \
  SvtxDeadMap.h \
  SvtxDeadMapv1.h \
  PHG4SvtxDeadArea.h \
  PHG4SvtxThresholds.h \
  PHG4SvtxDigitizer.h \
  PHG4SvtxClusterizer.h \
  SvtxHit.h \
  SvtxHit_v1.h \
  SvtxHitMap.h \
  SvtxHitMap_v1.h \
  SvtxCluster.h \
  SvtxCluster_v1.h \
  SvtxClusterMap.h \
  SvtxClusterMap_v1.h \
  SvtxTrackState.h \
  SvtxTrackState_v1.h \
  SvtxTrack.h \
  SvtxTrack_v1.h \
  SvtxTrack_FastSim.h \
  SvtxTrackMap.h \
  SvtxTrackMap_v1.h \
  SvtxBeamSpot.h \
  Cluster3D.h \
  Track3D.h \
  HelixTrackState.h \
  HelixHoughSpace.h \
  HelixHoughSpace_v1.h \
  HelixHoughBin.h \
  HelixHoughBin_v1.h \
  HelixHoughFuncs.h \
  HelixHoughFuncs_v1.h \
  HelixKalmanFilter.h \
  CellularAutomaton.h \
  CellularAutomaton_v1.h \
  FunctionGradHessian.h \
  SquareGradient.h \
  NewtonMinimizerGradHessian.h \
  VertexFitter.h \
  VertexFitFuncs.h \
  PHG4HoughTransform.h \
  PHG4HoughTransformTPC.h \
  PHG4KalmanPatRec.h \
  PHG4InitZVertexing.h \
  PHG4PatternReco.h \
  PHG4TrackGhostRejection.h \
  PHG4TrackKalmanFitter.h \
  PHG4TruthPatRec.h \
  PHG4TrackFastSim.h \
  PHRaveVertexing.h \
  PHG4SvtxMomentumRecal.h \
  PHG4GenFitTrackProjection.h \
  PHG4SvtxTrackProjection.h \
  PHG4SvtxBeamSpotReco.h \
  PHG4TPCClusterizer.h \
  SimpleTrack.h

libg4hough_io_la_SOURCES = \
  SvtxHit.C \
  SvtxHit_Dict.C \
  SvtxHit_v1.C \
  SvtxHit_v1_Dict.C \
  SvtxDeadMap.C \
  SvtxDeadMap_Dict.C \
  SvtxDeadMapv1.C \
  SvtxDeadMapv1_Dict.C \
  SvtxHitMap.C \
  SvtxHitMap_Dict.C \
  SvtxHitMap_v1.C \
  SvtxHitMap_v1_Dict.C \
  SvtxCluster.C \
  SvtxCluster_Dict.C \
  SvtxCluster_v1.C \
  SvtxCluster_v1_Dict.C \
  SvtxClusterMap.C \
  SvtxClusterMap_Dict.C \
  SvtxClusterMap_v1.C \
  SvtxClusterMap_v1_Dict.C \
  SvtxTrackState.C \
  SvtxTrackState_Dict.C \
  SvtxTrackState_v1.C \
  SvtxTrackState_v1_Dict.C \
  SvtxTrack.C \
  SvtxTrack_Dict.C \
  SvtxTrack_v1.C \
  SvtxTrack_FastSim.C \
  SvtxTrack_v1_Dict.C \
  SvtxTrack_FastSim_Dict.C \
  SvtxBeamSpot.C \
  SvtxBeamSpot_Dict.C \
  SvtxTrackMap.C \
  SvtxTrackMap_Dict.C \
  SvtxTrackMap_v1.C \
  SvtxTrackMap_v1_Dict.C \
  SvtxVertex.C \
  SvtxVertex_Dict.C \
  SvtxVertex_v1.C \
  SvtxVertex_v1_Dict.C \
  SvtxVertexMap.C \
  SvtxVertexMap_Dict.C \
  SvtxVertexMap_v1.C \
  SvtxVertexMap_v1_Dict.C \
  Cluster3D.C \
  Track3D.C \
  HelixTrackState.C \
  HelixHoughSpace.C \
  HelixHoughSpace_Dict.C \
  HelixHoughSpace_v1.C \
  HelixHoughSpace_v1_Dict.C \
  HelixHoughBin.C \
  HelixHoughBin_Dict.C \
  HelixHoughBin_v1.C \
  HelixHoughBin_v1_Dict.C \
  HelixHoughFuncs.C \
  HelixHoughFuncs_Dict.C \
  HelixHoughFuncs_v1.C \
  HelixHoughFuncs_v1_Dict.C \
  HelixKalmanFilter.C \
  CellularAutomaton.C \
  CellularAutomaton_Dict.C \
  CellularAutomaton_v1.C \
  CellularAutomaton_v1_Dict.C \
  SquareGradient.C \
  NewtonMinimizerGradHessian.C \
  VertexFitter.C \
  VertexFitFuncs.C

libg4hough_la_SOURCES = \
  PHG4SiliconTrackerDigitizer.C \
  PHG4SiliconTrackerDigitizer_Dict.C \
  PHG4SvtxDeadArea.C \
  PHG4SvtxDeadArea_Dict.C \
  PHG4SvtxDeadMapLoader.C \
  PHG4SvtxDeadMapLoader_Dict.C \
  PHG4SvtxThresholds.C \
  PHG4SvtxThresholds_Dict.C \
  PHG4SvtxDigitizer.C \
  PHG4SvtxDigitizer_Dict.C \
  PHG4SvtxClusterizer.C \
  PHG4SvtxClusterizer_Dict.C \
  PHG4TPCClusterizer.C \
  PHG4TPCClusterizer_Dict.C \
  PHG4HoughTransform.C \
  PHG4HoughTransform_Dict.C \
  PHG4HoughTransformTPC.C \
  PHG4HoughTransformTPC_Dict.C \
  PHG4KalmanPatRec.C \
  PHG4KalmanPatRec_Dict.C \
  PHRaveVertexing.C \
  PHRaveVertexing_Dict.C \
  PHG4InitZVertexing.C \
  PHG4InitZVertexing_Dict.C \
  PHG4PatternReco.C \
  PHG4PatternReco_Dict.C \
  PHG4TrackGhostRejection.C \
  PHG4TrackGhostRejection_Dict.C \
  PHG4TrackKalmanFitter.C \
  PHG4TruthPatRec.C \
  PHG4TrackFastSim.C \
  PHG4TrackKalmanFitter_Dict.C \
  PHG4TruthPatRec_Dict.C \
  PHG4TrackFastSim_Dict.C \
  PHG4SvtxMomentumRecal.C \
  PHG4SvtxMomentumRecal_Dict.C \
  PHG4GenFitTrackProjection.C \
  PHG4GenFitTrackProjection_Dict.C \
  PHG4SvtxTrackProjection.C \
  PHG4SvtxTrackProjection_Dict.C \
  PHG4SvtxBeamSpotReco.C \
  PHG4SvtxBeamSpotReco_Dict.C \
  SimpleTrack_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4hough_io \
  testexternals_g4hough

testexternals_g4hough_io_SOURCES = testexternals.C
testexternals_g4hough_io_LDADD = libg4hough_io.la

testexternals_g4hough_SOURCES = testexternals.C
testexternals_g4hough_LDADD = libg4hough.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* *.pcm
