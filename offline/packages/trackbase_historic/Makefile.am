AUTOMAKE_OPTIONS = foreign

# list of shared libraries to produce
lib_LTLIBRARIES = \
  libtrackbase_historic_io.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(ROOTSYS)/include \
  -I$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib

# sources for io library
libtrackbase_historic_io_la_SOURCES = \
  SvtxCluster.cc \
  SvtxCluster_Dict.cc \
  SvtxCluster_v1.cc \
  SvtxCluster_v1_Dict.cc \
  SvtxClusterMap.cc \
  SvtxClusterMap_Dict.cc \
  SvtxClusterMap_v1.cc \
  SvtxClusterMap_v1_Dict.cc \
  SvtxHit.cc \
  SvtxHit_Dict.cc \
  SvtxHit_v1.cc \
  SvtxHit_v1_Dict.cc \
  SvtxHitMap.cc \
  SvtxHitMap_Dict.cc \
  SvtxHitMap_v1.cc \
  SvtxHitMap_v1_Dict.cc \
  SvtxTrackState.cc \
  SvtxTrackState_Dict.cc \
  SvtxTrackState_v1.cc \
  SvtxTrackState_v1_Dict.cc \
  SvtxTrack.cc \
  SvtxTrack_Dict.cc \
  SvtxTrack_v1.cc \
  SvtxTrack_FastSim.cc \
  SvtxTrack_v1_Dict.cc \
  SvtxTrack_FastSim_Dict.cc \
  SvtxTrackMap.cc \
  SvtxTrackMap_Dict.cc \
  SvtxTrackMap_v1.cc \
  SvtxTrackMap_v1_Dict.cc \
  SvtxVertex.cc \
  SvtxVertex_Dict.cc \
  SvtxVertex_v1.cc \
  SvtxVertex_v1_Dict.cc \
  SvtxVertexMap.cc \
  SvtxVertexMap_Dict.cc \
  SvtxVertexMap_v1.cc \
  SvtxVertexMap_v1_Dict.cc 

libtrackbase_historic_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtrackbase_historic_io_la_LIBADD = \
  -lphool \
  -lCLHEP 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  SvtxHit.h \
  SvtxHitLinkDef.h \
  SvtxHit_v1.h \
  SvtxHit_v1LinkDef.h \
  SvtxHitMap.h \
  SvtxHitMapLinkDef.h \
  SvtxHitMap_v1.h \
  SvtxHitMap_v1LinkDef.h \
  SvtxCluster.h \
  SvtxClusterLinkDef.h \
  SvtxCluster_v1.h \
  SvtxCluster_v1LinkDef.h \
  SvtxClusterMap.h \
  SvtxClusterMapLinkDef.h \
  SvtxClusterMap_v1.h \
  SvtxClusterMap_v1LinkDef.h \
  SvtxTrack.h \
  SvtxTrack_v1.h \
  SvtxTrack_FastSim.h \
  SvtxTrackMap.h \
  SvtxTrackMap_v1.h \
  SvtxTrackState.h \
  SvtxTrackState_v1.h \
  SvtxVertex.h \
  SvtxVertex_v1.h \
  SvtxVertexMap.h \
  SvtxVertexMap_v1.h 


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
