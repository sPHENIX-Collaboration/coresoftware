AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcemc_io.la \
  libcemc.la

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4

RINCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include \
  -I${G4_MAIN}/include/Geant4

libcemc_io_la_SOURCES = \
  RawCluster_Dict.cc \
  RawClusterv1.cc \
  RawClusterv1_Dict.cc \
  RawClusterContainer.cc \
  RawClusterContainer_Dict.cc \
  RawTowerDefs_Dict.cc \
  RawTower_Dict.cc \
  RawTowerv1.cc \
  RawTowerv1_Dict.cc \
  RawTowerContainer.cc \
  RawTowerContainer_Dict.cc \
  RawTowerGeom.cc \
  RawTowerGeom_Dict.cc \
  RawTowerGeomv1.cc \
  RawTowerGeomv1_Dict.cc \
  RawTowerGeomv2.cc \
  RawTowerGeomv2_Dict.cc \
  RawTowerGeomv3.cc \
  RawTowerGeomv3_Dict.cc \
  RawTowerGeomContainer.cc \
  RawTowerGeomContainer_Dict.cc \
  RawTowerGeomContainerv1.cc \
  RawTowerGeomContainerv1_Dict.cc \
  RawTowerGeomContainer_Cylinderv1.cc \
  RawTowerGeomContainer_Cylinderv1_Dict.cc

libcemc_la_SOURCES = \
  Prototype2RawTowerBuilder.cc \
  Prototype2RawTowerBuilder_Dict.cc \
  RawClusterBuilder.cc \
  RawClusterBuilder_Dict.cc \
  RawClusterBuilderv1.cc \
  RawClusterBuilderv1_Dict.cc \
  RawClusterBuilderFwd.cc \
  RawClusterBuilderFwd_Dict.cc \
  RawTowerBuilder.cc \
  RawTowerBuilder_Dict.cc \
  RawTowerCombiner.cc \
  RawTowerCombiner_Dict.cc \
  RawTowerBuilderByHitIndex.cc \
  RawTowerBuilderByHitIndex_Dict.cc \
  RawTowerCalibration.cc \
  RawTowerCalibration_Dict.cc \
  RawTowerDigitizer.cc \
  RawTowerDigitizer_Dict.cc \
  BEmcCluster.cc \
  BEmcRec.cc

libcemc_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcemc_io_la_LIBADD = \
  -lphool

libcemc_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcemc_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4detectors \
  -lgsl \
  -lgslcblas \
  libcemc_io.la

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  RawCluster.h \
  RawClusterDefs.h \
  RawClusterContainer.h \
  RawTower.h \
  RawTowerDefs.h \
  RawTowerv1.h \
  RawTowerContainer.h  \
  RawTowerGeom.h \
  RawTowerGeomv1.h \
  RawTowerGeomv2.h \
  RawTowerGeomv3.h \
  RawTowerGeomContainer.h \
  RawTowerGeomContainerv1.h \
  RawTowerGeomContainer_Cylinderv1.h

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_cemc 

testexternals_cemc_SOURCES = testexternals.C
testexternals_cemc_LDADD = libcemc_io.la libcemc.la


testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
