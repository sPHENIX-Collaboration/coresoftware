AUTOMAKE_OPTIONS = foreign
# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtpctrackreco.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -fopenmp
AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  Tpc_ModuleTrackReco.h \
  Tpc_ModuleTrack.h \
  Tpc_ModuleTrackv1.h \
  Tpc_ModuleTrackContainer.h \
  Tpc_ModuleTrackContainerv1.h \
  Tpc_AssembledTrack.h \
  Tpc_AssembledTrackv1.h \
  Tpc_AssembledTrackContainer.h \
  Tpc_AssembledTrackContainerv1.h \
  Tpc_AssembledTrackReco.h \
  Tpc_PolyTrack.h \
  Tpc_PolyTrackv1.h \
  Tpc_PolyTrackContainer.h \
  Tpc_PolyTrackContainerv1.h \
  Tpc_PolyTrackVertex.h \
  Tpc_PolyTrackVertexv1.h \
  Tpc_PolyTrackVertexContainer.h \
  Tpc_PolyTrackVertexContainerv1.h \
  Tpc_PolyTrackVertexer.h \
  Tpc_PolyCluster.h \
  Tpc_PolyClusterv1.h \
  Tpc_PolyClusterContainer.h \
  Tpc_PolyClusterContainerv1.h \
  Tpc_PolyClusterizer.h \
  Tpc_PolyTrackReco.h \
  Tpc_FittingTools.h \
  IdealPadMap.h


ROOTDICTS = \
  Tpc_AssembledTrack_Dict.cc \
  Tpc_AssembledTrackContainer_Dict.cc \
  Tpc_AssembledTrackContainerv1_Dict.cc \
  Tpc_AssembledTrackv1_Dict.cc \
  Tpc_ModuleTrack_Dict.cc \
  Tpc_ModuleTrackContainer_Dict.cc \
  Tpc_ModuleTrackContainerv1_Dict.cc \
  Tpc_ModuleTrackv1_Dict.cc \
  Tpc_PolyCluster_Dict.cc \
  Tpc_PolyClusterContainer_Dict.cc \
  Tpc_PolyClusterContainerv1_Dict.cc \
  Tpc_PolyClusterv1_Dict.cc \
  Tpc_PolyTrack_Dict.cc \
  Tpc_PolyTrackContainer_Dict.cc \
  Tpc_PolyTrackContainerv1_Dict.cc \
  Tpc_PolyTrackVertex_Dict.cc \
  Tpc_PolyTrackVertexContainer_Dict.cc \
  Tpc_PolyTrackVertexContainerv1_Dict.cc \
  Tpc_PolyTrackVertexv1_Dict.cc \
  Tpc_PolyTrackv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  Tpc_AssembledTrack_Dict_rdict.pcm \
  Tpc_AssembledTrackContainer_Dict_rdict.pcm \
  Tpc_AssembledTrackContainerv1_Dict_rdict.pcm \
  Tpc_AssembledTrackv1_Dict_rdict.pcm \
  Tpc_ModuleTrack_Dict_rdict.pcm \
  Tpc_ModuleTrackContainer_Dict_rdict.pcm \
  Tpc_ModuleTrackContainerv1_Dict_rdict.pcm \
  Tpc_ModuleTrackv1_Dict_rdict.pcm \
  Tpc_PolyCluster_Dict_rdict.pcm \
  Tpc_PolyClusterContainer_Dict_rdict.pcm \
  Tpc_PolyClusterContainerv1_Dict_rdict.pcm \
  Tpc_PolyClusterv1_Dict_rdict.pcm \
  Tpc_PolyTrack_Dict_rdict.pcm \
  Tpc_PolyTrackContainer_Dict_rdict.pcm \
  Tpc_PolyTrackContainerv1_Dict_rdict.pcm \
  Tpc_PolyTrackVertex_Dict_rdict.pcm \
  Tpc_PolyTrackVertexContainer_Dict_rdict.pcm \
  Tpc_PolyTrackVertexContainerv1_Dict_rdict.pcm \
  Tpc_PolyTrackVertexv1_Dict_rdict.pcm \
  Tpc_PolyTrackv1_Dict_rdict.pcm

libtpctrackreco_la_SOURCES = \
  Tpc_ModuleTrackReco.cc \
  Tpc_ModuleTrackv1.cc \
  Tpc_ModuleTrackContainerv1.cc \
  Tpc_AssembledTrackv1.cc \
  Tpc_AssembledTrackContainerv1.cc \
  Tpc_AssembledTrackReco.cc \
  Tpc_PolyTrackv1.cc \
  Tpc_PolyTrackContainerv1.cc \
  Tpc_PolyTrackVertexv1.cc \
  Tpc_PolyTrackVertexContainerv1.cc \
  Tpc_PolyTrackVertexer.cc \
  Tpc_PolyClusterv1.cc \
  Tpc_PolyClusterContainerv1.cc \
  Tpc_PolyClusterizer.cc \
  Tpc_PolyTrackReco.cc \
  Tpc_FittingTools.cc \
  IdealPadMap.cc \
  $(ROOTDICTS)

libtpctrackreco_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lffamodules \
  -lcdbobjects \
  -lg4eval \
  -lg4eval_io \
  -lphg4hit \
  -lmbd_io \
  -lg4detectors_io \
  -lmvtx \
  -lmvtx_io \
  -lintt \
  -lintt_io \
  -lmicromegas_io \
  -lSubsysReco \
  -ltpc_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
  -lphfield \
  -lPHGarfield

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_tpc_track_reco


testexternals_tpc_track_reco_SOURCES = testexternals.cc
testexternals_tpc_track_reco_LDADD = libtpctrackreco.la

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