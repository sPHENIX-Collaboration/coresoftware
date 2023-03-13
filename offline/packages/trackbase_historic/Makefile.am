##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrackbase_historic_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  ActsTransformations.h \
  TrackSeed.h \
  TrackSeed_v1.h \
  SvtxTrackSeed_v1.h \
  TrackSeed_FastSim_v1.h \
  TrackSeedContainer.h \
  TrackSeedContainer_v1.h \
  PHG4ParticleSvtxMap.h \
  PHG4ParticleSvtxMap_v1.h \
  SvtxPHG4ParticleMap.h \
  SvtxPHG4ParticleMap_v1.h \
  SvtxTrack.h \
  SvtxTrack_v1.h \
  SvtxTrack_v2.h \
  SvtxTrack_v3.h \
  SvtxTrack_v4.h \
  SvtxTrack_FastSim.h \
  SvtxTrack_FastSim_v1.h \
  SvtxTrack_FastSim_v2.h \
  SvtxTrack_FastSim_v3.h \
  SvtxAlignmentStateMap.h \
  SvtxAlignmentStateMap_v1.h \
  SvtxTrackMap.h \
  SvtxTrackMap_v1.h \
  SvtxTrackMap_v2.h \
  SvtxTrackCaloClusterMap.h \
  SvtxTrackCaloClusterMap_v1.h \
  SvtxAlignmentState.h \
  SvtxAlignmentState_v1.h \
  SvtxTrackState.h \
  SvtxTrackState_v1.h \
  SvtxVertex.h \
  SvtxVertex_v1.h \
  SvtxVertexMap.h \
  SvtxVertexMap_v1.h 

ROOTDICTS = \
  TrackSeed_Dict.cc \
  TrackSeed_v1_Dict.cc \
  SvtxTrackSeed_v1_Dict.cc \
  TrackSeed_FastSim_v1_Dict.cc \
  TrackSeedContainer_Dict.cc \
  TrackSeedContainer_v1_Dict.cc \
  PHG4ParticleSvtxMap_Dict.cc \
  PHG4ParticleSvtxMap_v1_Dict.cc \
  SvtxPHG4ParticleMap_Dict.cc \
  SvtxPHG4ParticleMap_v1_Dict.cc \
  SvtxAlignmentState_Dict.cc \
  SvtxAlignmentState_v1_Dict.cc \
  SvtxTrack_Dict.cc \
  SvtxTrackState_Dict.cc \
  SvtxTrackState_v1_Dict.cc \
  SvtxTrack_v1_Dict.cc \
  SvtxTrack_v2_Dict.cc \
  SvtxTrack_v3_Dict.cc \
  SvtxTrack_v4_Dict.cc \
  SvtxTrack_FastSim_Dict.cc \
  SvtxTrack_FastSim_v1_Dict.cc \
  SvtxTrack_FastSim_v2_Dict.cc \
  SvtxTrack_FastSim_v3_Dict.cc \
  SvtxAlignmentStateMap_Dict.cc \
  SvtxAlignmentStateMap_v1_Dict.cc \
  SvtxTrackMap_Dict.cc \
  SvtxTrackMap_v1_Dict.cc \
  SvtxTrackMap_v2_Dict.cc \
  SvtxTrackCaloClusterMap_Dict.cc \
  SvtxTrackCaloClusterMap_v1_Dict.cc \
  SvtxVertex_Dict.cc \
  SvtxVertex_v1_Dict.cc \
  SvtxVertexMap_Dict.cc \
  SvtxVertexMap_v1_Dict.cc 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TrackSeed_Dict_rdict.pcm \
  TrackSeed_v1_Dict_rdict.pcm \
  SvtxTrackSeed_v1_Dict_rdict.pcm \
  TrackSeed_FastSim_v1_Dict_rdict.pcm \
  TrackSeedContainer_Dict_rdict.pcm \
  TrackSeedContainer_v1_Dict_rdict.pcm \
  PHG4ParticleSvtxMap_Dict_rdict.pcm \
  PHG4ParticleSvtxMap_v1_Dict_rdict.pcm \
  SvtxPHG4ParticleMap_Dict_rdict.pcm \
  SvtxPHG4ParticleMap_v1_Dict_rdict.pcm \
  SvtxAlignmentState_Dict_rdict.pcm \
  SvtxAlignmentState_v1_Dict_rdict.pcm \
  SvtxTrack_Dict_rdict.pcm \
  SvtxTrackState_Dict_rdict.pcm \
  SvtxTrackState_v1_Dict_rdict.pcm \
  SvtxTrack_v1_Dict_rdict.pcm \
  SvtxTrack_v2_Dict_rdict.pcm \
  SvtxTrack_v3_Dict_rdict.pcm \
  SvtxTrack_v4_Dict_rdict.pcm \
  SvtxTrack_FastSim_Dict_rdict.pcm \
  SvtxTrack_FastSim_v1_Dict_rdict.pcm \
  SvtxTrack_FastSim_v2_Dict_rdict.pcm \
  SvtxTrack_FastSim_v3_Dict_rdict.pcm \
  SvtxAlignmentStateMap_Dict_rdict.pcm \
  SvtxAlignmentStateMap_v1_Dict_rdict.pcm \
  SvtxTrackMap_Dict_rdict.pcm \
  SvtxTrackMap_v1_Dict_rdict.pcm \
  SvtxTrackMap_v2_Dict_rdict.pcm \
  SvtxTrackCaloClusterMap_Dict_rdict.pcm \
  SvtxTrackCaloClusterMap_v1_Dict_rdict.pcm \
  SvtxVertex_Dict_rdict.pcm \
  SvtxVertex_v1_Dict_rdict.pcm \
  SvtxVertexMap_Dict_rdict.pcm \
  SvtxVertexMap_v1_Dict_rdict.pcm 

# sources for io library
libtrackbase_historic_io_la_SOURCES = \
  $(ROOTDICTS) \
  ActsTransformations.cc \
  TrackSeed.cc \
  TrackSeed_v1.cc \
  SvtxTrackSeed_v1.cc \
  TrackSeed_FastSim_v1.cc \
  TrackSeedContainer.cc \
  TrackSeedContainer_v1.cc \
  PHG4ParticleSvtxMap.cc \
  PHG4ParticleSvtxMap_v1.cc \
  SvtxPHG4ParticleMap.cc \
  SvtxPHG4ParticleMap_v1.cc \
  SvtxAlignmentState.cc \
  SvtxAlignmentState_v1.cc \
  SvtxTrackState_v1.cc \
  SvtxTrack.cc \
  SvtxTrack_v1.cc \
  SvtxTrack_v2.cc \
  SvtxTrack_v3.cc \
  SvtxTrack_v4.cc \
  SvtxTrack_FastSim.cc \
  SvtxTrack_FastSim_v1.cc \
  SvtxTrack_FastSim_v2.cc \
  SvtxTrack_FastSim_v3.cc \
  SvtxAlignmentStateMap.cc \
  SvtxAlignmentStateMap_v1.cc \
  SvtxTrackMap.cc \
  SvtxTrackMap_v1.cc \
  SvtxTrackMap_v2.cc \
  SvtxTrackCaloClusterMap.cc \
  SvtxTrackCaloClusterMap_v1.cc \
  SvtxVertex.cc \
  SvtxVertex_v1.cc \
  SvtxVertexMap.cc \
  SvtxVertexMap_v1.cc

libtrackbase_historic_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libtrackbase_historic_io_la_LIBADD = \
  -lphool \
  -lActsCore \
  -ltrack_io

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;



################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_trackbase_historic_io

testexternals_trackbase_historic_io_SOURCES = testexternals.cc
testexternals_trackbase_historic_io_LDADD = libtrackbase_historic_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
