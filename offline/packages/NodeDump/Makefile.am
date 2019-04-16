AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libphnodedump.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

pkginclude_HEADERS = \
  Dumper.h \
  DumpObject.h \
  PHNodeDump.h

if ! MAKEROOT6
 ROOT5_DICTS = \
  Dumper_Dict.cc \
  PHNodeDump_Dict.cc
endif

libphnodedump_la_SOURCES = \
  $(ROOT5_DICTS) \
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
  DumpTowerBackground.cc \
  DumpVariableArray.cc \
  DumpObject.cc \
  PHNodeDump.cc

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libphnodedump_la_LIBADD = \
  -lg4bbc_io \
  -lcalo_io \
  -lcalotrigger_io \
  -lffaobjects \
  -lg4detectors_io \
  -lg4hough_io \
  -lg4jets_io \
  -lphg4hit \
  -lSubsysReco \
  -ltrackbase_historic_io \
  -lvararray


# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

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
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
