##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcalo_io.la

libcalo_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalo_io_la_LIBADD = \
  -lphool

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  RawClusterUtility.h \
  RawCluster.h \
  RawClusterv1.h \
  RawClusterDefs.h \
  RawClusterContainer.h \
  RawTower.h \
  RawTowerDefs.h \
  RawTowerv1.h \
  RawTowerv2.h \
  RawTowerContainer.h  \
  RawTowerDeadMap.h  \
  RawTowerDeadMapv1.h  \
  RawTowerGeom.h \
  RawTowerGeomv1.h \
  RawTowerGeomv2.h \
  RawTowerGeomv3.h \
  RawTowerGeomv4.h \
  RawTowerGeomContainer.h \
  RawTowerGeomContainerv1.h \
  RawTowerGeomContainer_Cylinderv1.h \
  TowerInfo.h \
  TowerInfov1.h \
  TowerInfoContainer.h \
  TowerInfoContainerv1.h

ROOTDICTS = \
  RawCluster_Dict.cc \
  RawClusterv1_Dict.cc \
  RawClusterContainer_Dict.cc \
  RawTower_Dict.cc \
  RawTowerv1_Dict.cc \
  RawTowerv2_Dict.cc \
  RawTowerContainer_Dict.cc \
  RawTowerDeadMap_Dict.cc \
  RawTowerGeom_Dict.cc \
  RawTowerDeadMapv1_Dict.cc \
  RawTowerGeomv1_Dict.cc \
  RawTowerGeomv2_Dict.cc \
  RawTowerGeomv3_Dict.cc \
  RawTowerGeomv4_Dict.cc \
  RawTowerGeomContainer_Dict.cc \
  RawTowerGeomContainerv1_Dict.cc \
  RawTowerGeomContainer_Cylinderv1_Dict.cc \
  TowerInfo_Dict.cc \
  TowerInfov1_Dict.cc \
  TowerInfoContainer_Dict.cc \
  TowerInfoContainerv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  RawCluster_Dict_rdict.pcm \
  RawClusterv1_Dict_rdict.pcm \
  RawClusterContainer_Dict_rdict.pcm \
  RawTower_Dict_rdict.pcm \
  RawTowerv1_Dict_rdict.pcm \
  RawTowerv2_Dict_rdict.pcm \
  RawTowerContainer_Dict_rdict.pcm \
  RawTowerDeadMap_Dict_rdict.pcm \
  RawTowerGeom_Dict_rdict.pcm \
  RawTowerDeadMapv1_Dict_rdict.pcm \
  RawTowerGeomv1_Dict_rdict.pcm \
  RawTowerGeomv2_Dict_rdict.pcm \
  RawTowerGeomv3_Dict_rdict.pcm \
  RawTowerGeomv4_Dict_rdict.pcm \
  RawTowerGeomContainer_Dict_rdict.pcm \
  RawTowerGeomContainerv1_Dict_rdict.pcm \
  RawTowerGeomContainer_Cylinderv1_Dict_rdict.pcm \
  TowerInfo_Dict_rdict.pcm \
  TowerInfov1_Dict_rdict.pcm \
  TowerInfoContainer_Dict_rdict.pcm \
  TowerInfoContainerv1_Dict_rdict.pcm

libcalo_io_la_SOURCES = \
  $(ROOTDICTS) \
  RawCluster.cc \
  RawClusterv1.cc \
  RawClusterContainer.cc \
  RawTower.cc \
  RawTowerv1.cc \
  RawTowerv2.cc \
  RawTowerContainer.cc \
  RawTowerDeadMap.cc \
  RawTowerDeadMapv1.cc \
  RawTowerGeom.cc \
  RawTowerGeomv1.cc \
  RawTowerGeomv2.cc \
  RawTowerGeomv3.cc \
  RawTowerGeomv4.cc \
  RawTowerGeomContainer.cc \
  RawTowerGeomContainerv1.cc \
  RawTowerGeomContainer_Cylinderv1.cc \
  TowerInfov1.cc \
  TowerInfoContainer.cc \
  TowerInfoContainerv1.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_calo_io

BUILT_SOURCES = testexternals.cc

testexternals_calo_io_SOURCES = testexternals.cc
testexternals_calo_io_LDADD = libcalo_io.la

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
