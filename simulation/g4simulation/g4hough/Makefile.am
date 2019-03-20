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
 -Werror -DRAVE -DRaveDllExport=

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libg4hough_io_la_LIBADD = \
  -lphool \
  -lg4detectors_io \
  -lintt \
  -lmvtx

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
  -lg4intt \
  -lg4mvtx \
  -lmvtx \
  -lintt \
  -ltrackbase_historic_io \
  libg4hough_io.la

pkginclude_HEADERS = \
  CellularAutomaton.h \
  CellularAutomaton_v1.h \
  Cluster3D.h \
  FunctionGradHessian.h \
  HelixHoughSpace.h \
  HelixHoughSpace_v1.h \
  HelixHoughBin.h \
  HelixHoughBin_v1.h \
  HelixHoughFuncs.h \
  HelixHoughFuncs_v1.h \
  HelixKalmanFilter.h \
  HelixTrackState.h \
  NewtonMinimizerGradHessian.h \
  PHG4GenFitTrackProjection.h \
  PHG4HoughTransform.h \
  PHG4HoughTransformTPC.h \
  PHG4InitZVertexing.h \
  PHG4KalmanPatRec.h \
  PHG4PatternReco.h \
  PHG4SvtxBeamSpotReco.h \
  PHG4SvtxDeadArea.h \
  PHG4SvtxMomentumRecal.h \
  PHG4SvtxThresholds.h \
  PHG4SvtxTrackProjection.h \
  PHG4TrackGhostRejection.h \
  PHG4TrackKalmanFitter.h \
  PHG4TruthPatRec.h \
  PHRaveVertexing.h \
  SimpleTrack.h \
  SquareGradient.h \
  SvtxBeamSpot.h \
  Track3D.h \
  VertexFitter.h \
  VertexFitFuncs.h

ROOTDICTS = \
  HelixHoughSpace_Dict.cc \
  HelixHoughSpace_v1_Dict.cc \
  HelixHoughBin_Dict.cc \
  HelixHoughBin_v1_Dict.cc \
  HelixHoughFuncs_Dict.cc \
  HelixHoughFuncs_v1_Dict.cc \
  CellularAutomaton_Dict.cc \
  CellularAutomaton_v1_Dict.cc \
  SimpleTrack_Dict.cc \
  SvtxBeamSpot_Dict.cc

if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  HelixHoughSpace_Dict_rdict.pcm \
  HelixHoughSpace_v1_Dict_rdict.pcm \
  HelixHoughBin_Dict_rdict.pcm \
  HelixHoughBin_v1_Dict_rdict.pcm \
  HelixHoughFuncs_Dict_rdict.pcm \
  HelixHoughFuncs_v1_Dict_rdict.pcm \
  CellularAutomaton_Dict_rdict.pcm \
  CellularAutomaton_v1_Dict_rdict.pcm \
  SimpleTrack_Dict_rdict.pcm \
  SvtxBeamSpot_Dict_rdict.pcm
else
  ROOT5_DICTS = \
  PHG4SvtxDeadArea_Dict.cc \
  PHG4SvtxThresholds_Dict.cc \
  PHG4HoughTransform_Dict.cc \
  PHG4HoughTransformTPC_Dict.cc \
  PHG4KalmanPatRec_Dict.cc \
  PHRaveVertexing_Dict.cc \
  PHG4InitZVertexing_Dict.cc \
  PHG4TruthPatRec_Dict.cc \
  PHG4PatternReco_Dict.cc \
  PHG4TrackGhostRejection_Dict.cc \
  PHG4TrackKalmanFitter_Dict.cc \
  PHG4SvtxMomentumRecal_Dict.cc \
  PHG4GenFitTrackProjection_Dict.cc \
  PHG4SvtxTrackProjection_Dict.cc \
  PHG4SvtxBeamSpotReco_Dict.cc
endif

libg4hough_io_la_SOURCES = \
  $(ROOTDICTS) \
  Cluster3D.C \
  Track3D.C \
  HelixTrackState.C \
  HelixHoughSpace_v1.C \
  HelixHoughBin_v1.C \
  HelixHoughFuncs_v1.C \
  HelixKalmanFilter.C \
  CellularAutomaton_v1.C \
  SquareGradient.C \
  SvtxBeamSpot.C \
  NewtonMinimizerGradHessian.C \
  VertexFitter.C \
  VertexFitFuncs.C

libg4hough_la_SOURCES = \
  $(ROOT5_DICTS) \
  PHG4SvtxDeadArea.C \
  PHG4SvtxThresholds.C \
  PHG4HoughTransform.C \
  PHG4HoughTransformTPC.C \
  PHG4KalmanPatRec.C \
  PHRaveVertexing.C \
  PHG4InitZVertexing.C \
  PHG4PatternReco.C \
  PHG4TrackGhostRejection.C \
  PHG4TrackKalmanFitter.C \
  PHG4TruthPatRec.C \
  PHG4SvtxMomentumRecal.C \
  PHG4GenFitTrackProjection.C \
  PHG4SvtxTrackProjection.C \
  PHG4SvtxBeamSpotReco.C

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4hough_io \
  testexternals_g4hough

BUILT_SOURCES = testexternals.cc

testexternals_g4hough_io_SOURCES = testexternals.cc
testexternals_g4hough_io_LDADD = libg4hough_io.la

testexternals_g4hough_SOURCES = testexternals.cc
testexternals_g4hough_LDADD = libg4hough.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
