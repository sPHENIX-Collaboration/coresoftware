AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.cc

lib_LTLIBRARIES = libphnodedump.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

pkginclude_HEADERS = \
  Dumper.h \
  DumpObject.h \
  PHNodeDump.h

libphnodedump_la_SOURCES = \
  Dumper.cc \
  DumpBbcVertexMap.cc \
  DumpCaloTriggerInfo.cc \
  DumpGlobalVertexMap.cc \
  DumpJetMap.cc \
  DumpPdbParameterMap.cc \
  DumpPdbParameterMapContainer.cc \
  DumpPHG4BlockGeomContainer.cc \
  DumpPHG4BlockCellGeomContainer.cc \
  DumpPHG4CellContainer.cc \
  DumpPHG4CylinderCellContainer.cc \
  DumpPHG4CylinderCellGeomContainer.cc \
  DumpPHG4CylinderGeomContainer.cc \
  DumpPHG4HitContainer.cc \
  DumpPHG4ScintillatorSlatContainer.cc \
  DumpPHG4TruthInfoContainer.cc \
  DumpPHG4InEvent.cc \
  DumpRawClusterContainer.cc \
  DumpRawTowerContainer.cc \
  DumpRawTowerGeomContainer.cc \
  DumpRunHeader.cc \
  DumpSvtxClusterMap.cc \
  DumpSvtxHitMap.cc \
  DumpSvtxTrackMap.cc \
  DumpSvtxVertexMap.cc \
  DumpSyncObject.cc \
  DumpVariableArray.cc \
  DumpObject.cc \
  PHNodeDump.cc \
  PHNodeDumpDict.cc

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libphnodedump_la_LIBADD = \
  -lg4bbc_io \
  -lcalo_util \
  -lcalotrigger_io \
  -lffaobjects \
  -lg4detectors \
  -lg4hough_io \
  -lg4jets_io \
  -lg4testbench \
  -lphg4hit \
  -lvararray

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libphnodedump.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

PHNodeDumpDict.cc : \
  PHNodeDump.h \
  Dumper.h \
  PHNodeDumpLinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*
