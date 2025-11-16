##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcalo_io.la

libcalo_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

# List of shared libraries to produce
if USE_ONLINE
pkginclude_HEADERS = \
  RawTowerDefs.h \
  TowerInfoDefs.h

libcalo_io_la_SOURCES = \
  TowerInfoDefs.cc

else

libcalo_io_la_LIBADD = \
  -lphool

pkginclude_HEADERS = \
  PhotonClusterv1.h \
  RawClusterUtility.h \
  RawCluster.h \
  RawClusterv1.h \
  RawClusterv2.h \
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
  RawTowerGeomv5.h \
  RawTowerGeomContainer.h \
  RawTowerGeomContainerv1.h \
  RawTowerGeomContainer_Cylinderv1.h \
  TowerInfoDefs.h \
  TowerInfo.h \
  TowerInfov1.h \
  TowerInfov2.h \
  TowerInfov3.h \
  TowerInfov4.h \
  TowerInfoSimv1.h \
  TowerInfoSimv2.h \
  TowerInfoContainer.h \
  TowerInfoContainerv1.h \
  TowerInfoContainerv2.h \
  TowerInfoContainerv3.h \
  TowerInfoContainerv4.h \
  TowerInfoContainerSimv1.h \
  TowerInfoContainerSimv2.h

ROOTDICTS = \
  PhotonClusterv1_Dict.cc \
  RawCluster_Dict.cc \
  RawClusterv1_Dict.cc \
  RawClusterv2_Dict.cc \
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
  RawTowerGeomv5_Dict.cc \
  RawTowerGeomContainer_Dict.cc \
  RawTowerGeomContainerv1_Dict.cc \
  RawTowerGeomContainer_Cylinderv1_Dict.cc \
  TowerInfo_Dict.cc \
  TowerInfov1_Dict.cc \
  TowerInfov2_Dict.cc \
  TowerInfov3_Dict.cc \
  TowerInfov4_Dict.cc \
  TowerInfoSimv1_Dict.cc \
  TowerInfoSimv2_Dict.cc \
  TowerInfoContainer_Dict.cc \
  TowerInfoContainerv1_Dict.cc \
  TowerInfoContainerv2_Dict.cc \
  TowerInfoContainerv3_Dict.cc \
  TowerInfoContainerv4_Dict.cc \
  TowerInfoContainerSimv1_Dict.cc \
  TowerInfoContainerSimv2_Dict.cc

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libcalo_io_la_SOURCES = \
  $(ROOTDICTS) \
  PhotonClusterv1.cc \
  RawCluster.cc \
  RawClusterv1.cc \
  RawClusterv2.cc \
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
  RawTowerGeomv5.cc \
  RawTowerGeomContainer.cc \
  RawTowerGeomContainerv1.cc \
  RawTowerGeomContainer_Cylinderv1.cc \
  TowerInfov1.cc \
  TowerInfov2.cc \
  TowerInfov3.cc \
  TowerInfov4.cc \
  TowerInfoSimv1.cc \
  TowerInfoSimv2.cc \
  TowerInfoDefs.cc \
  TowerInfoContainer.cc \
  TowerInfoContainerv1.cc \
  TowerInfoContainerv2.cc \
  TowerInfoContainerv3.cc \
  TowerInfoContainerv4.cc \
  TowerInfoContainerSimv1.cc \
  TowerInfoContainerSimv2.cc
endif

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
