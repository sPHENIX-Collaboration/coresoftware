AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libphnodedump.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  Dumper.h \
  DumpObject.h \
  PHNodeDump.h

libphnodedump_la_SOURCES = \
  Dumper.cc \
  DumpBbcVertexMap.cc \
  DumpCaloTriggerInfo.cc \
  DumpCdbUrlSave.cc \
  DumpCentralityInfo.cc \
  DumpEpInfo.cc \
  DumpEventHeader.cc \
  DumpFlagSave.cc \
  DumpGlobalVertexMap.cc \
  DumpInttDeadMap.cc \
  DumpJetMap.cc \
  DumpParticleFlowElementContainer.cc \
  DumpPdbParameterMap.cc \
  DumpPdbParameterMapContainer.cc \
  DumpPHFieldConfig.cc \
  DumpPHG4BlockGeomContainer.cc \
  DumpPHG4BlockCellGeomContainer.cc \
  DumpPHG4CellContainer.cc \
  DumpPHG4CylinderCellContainer.cc \
  DumpPHG4CylinderCellGeomContainer.cc \
  DumpPHG4CylinderGeomContainer.cc \
  DumpPHG4HitContainer.cc \
  DumpPHG4ParticleSvtxMap.cc \
  DumpPHG4ScintillatorSlatContainer.cc \
  DumpPHG4TruthInfoContainer.cc \
  DumpPHG4InEvent.cc \
  DumpPHGenIntegral.cc \
  DumpPHHepMCGenEventMap.cc \
  DumpRawClusterContainer.cc \
  DumpRawTowerContainer.cc \
  DumpRawTowerGeomContainer.cc \
  DumpRunHeader.cc \
  DumpSvtxPHG4ParticleMap.cc \
  DumpSvtxTrackMap.cc \
  DumpSvtxVertexMap.cc \
  DumpSyncObject.cc \
  DumpTowerBackground.cc \
  DumpTpcSeedTrackMap.cc \
  DumpTrackSeedContainer.cc \
  DumpTrkrClusterContainer.cc \
  DumpTrkrClusterHitAssoc.cc \
  DumpTrkrHitSetContainer.cc \
  DumpTrkrHitTruthClusters.cc \
  DumpTrkrClusterCrossingAssoc.cc \
  DumpTrkrHitTruthAssoc.cc \
  DumpVariableArray.cc \
  DumpObject.cc \
  PHNodeDump.cc

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libphnodedump_la_LIBADD = \
  -lg4bbc_io \
  -lcalo_io \
  -lcalotrigger_io \
  -leventplane_io \
  -lffaobjects \
  -lg4detectors_io \
  -lphfield_io \
  -lg4jets_io \
  -lHepMC \
  -lparticleflow_io \
  -lphg4hit \
  -lSubsysReco \
  -ltrackbase_historic_io \
  -ltrack_io \
  -lphhepmc_io\
  -lvararray


noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libphnodedump.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
