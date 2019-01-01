AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtrack_reco.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I${G4_MAIN}/include \
  -I$(ROOTSYS)/include \
  -I$(OFFLINE_MAIN)/include/eigen3

AM_CXXFLAGS = \
 -Werror -msse2 -DRAVE -DRaveDllExport=
   
libtrack_reco_la_SOURCES = \
  AssocInfoContainer.C \
  AssocInfoContainer_Dict.C \
  PHInitVertexing.C \
  PHInitVertexing_Dict.C \
  PHTrackSeeding.C \
  PHTrackSeeding_Dict.C \
  PHTrackPropagating.C \
  PHTrackPropagating_Dict.C \
  PHTrackFitting.C \
  PHTrackFitting_Dict.C \
  PH3DVertexing.C \
  PH3DVertexing_Dict.C \
  PHTruthVertexing.C \
  PHTruthVertexing_Dict.C \
  PHTruthTrackSeeding.C \
  PHTruthTrackSeeding_Dict.C \
  PHHoughSeeding.C \
  PHHoughSeeding_Dict.C \
  PHGenFitTrkProp.C \
  PHGenFitTrkProp_Dict.C \
  PHKalmanTrkFitter.C \
  PHKalmanTrkFitter_Dict.C \
  PHRaveVertexing.C \
  PHRaveVertexing_Dict.C

libtrack_reco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtrack_reco_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4testbench \
  -lg4detectors \
  -lgsl \
  -lgslcblas \
  -lHelixHough \
  -lFitNewton \
  -lSeamstress \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
  -lg4bbc_io \
  -ltrackbase_historic \
  -lphparameter
  
# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

# pkginclude_HEADERS = 

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_track_reco


testexternals_track_reco_SOURCES = testexternals.C
testexternals_track_reco_LDADD = libtrack_reco.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
