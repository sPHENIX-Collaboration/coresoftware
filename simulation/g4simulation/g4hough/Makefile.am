AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4 \
  -I`root-config --incdir` \
  -I$(OFFLINE_MAIN)/include/eigen3

lib_LTLIBRARIES = \
   libg4hough_io.la \
   libg4hough.la

AM_CXXFLAGS = \
 -Wall -Werror -msse2

libg4hough_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool

libg4hough_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `geant4-config --libs` \
  -lg4detectors \
  -lfun4all \
  -lHelixHough \
  -lFitNewton \
  -lSeamstress \
  -lcemc \
  -lg4hough_io

pkginclude_HEADERS = \
  SvtxVertex.h \
  SvtxVertexMap.h \
  PHG4SvtxDeadArea.h \
  PHG4SvtxThresholds.h \
  PHG4SvtxAddConnectedCells.h \
  PHG4SvtxDigitizer.h \
  PHG4SvtxClusterizer.h \
  SvtxHit.h \
  SvtxHitMap.h \
  SvtxCluster.h \
  SvtxClusterMap.h \
  SvtxTrack.h \
  SvtxTrackMap.h \
  SvtxBeamSpot.h \
  PHG4HoughTransform.h \
  PHG4TrackGhostRejection.h \
  PHG4SvtxTrackProjection.h \
  PHG4SvtxBeamSpotReco.h

#pkginclude_HEADERS = $(include_HEADERS)

libg4hough_io_la_SOURCES = \
  SvtxHit.C \
  SvtxHit_Dict.C \
  SvtxHitMap.C \
  SvtxHitMap_Dict.C \
  SvtxCluster.C \
  SvtxCluster_Dict.C \
  SvtxClusterMap.C \
  SvtxClusterMap_Dict.C \
  SvtxTrack.C \
  SvtxTrack_Dict.C \
  SvtxBeamSpot.C \
  SvtxBeamSpot_Dict.C \
  SvtxTrackMap.C \
  SvtxTrackMap_Dict.C \
  SvtxVertex.C \
  SvtxVertex_Dict.C \
  SvtxVertexMap.C \
  SvtxVertexMap_Dict.C

libg4hough_la_SOURCES = \
  PHG4SvtxDeadArea.C \
  PHG4SvtxDeadArea_Dict.C \
  PHG4SvtxThresholds.C \
  PHG4SvtxThresholds_Dict.C \
  PHG4SvtxAddConnectedCells.C \
  PHG4SvtxAddConnectedCells_Dict.C \
  PHG4SvtxDigitizer.C \
  PHG4SvtxDigitizer_Dict.C \
  PHG4SvtxClusterizer.C \
  PHG4SvtxClusterizer_Dict.C \
  PHG4HoughTransform.C \
  PHG4HoughTransform_Dict.C \
  PHG4TrackGhostRejection.C \
  PHG4TrackGhostRejection_Dict.C \
  PHG4SvtxTrackProjection.C \
  PHG4SvtxTrackProjection_Dict.C \
  PHG4SvtxBeamSpotReco.C \
  PHG4SvtxBeamSpotReco_Dict.C

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

%_Dict.cpp: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libg4hough_io.la \
  libg4hough.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict*
