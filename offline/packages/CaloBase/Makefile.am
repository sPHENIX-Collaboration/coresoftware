AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcalo_io.la \
  libcalo_util.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

libcalo_io_la_SOURCES = \
  RawCluster.cc \
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

libcalo_util_la_SOURCES = \
  RawClusterUtility.cc \
  RawClusterUtility_Dict.cc

libcalo_io_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalo_io_la_LIBADD = \
  -lphool 

libcalo_util_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalo_util_la_LIBADD = \
  -lCLHEP \
  libcalo_io.la

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  RawClusterUtility.h \
  RawCluster.h \
  RawClusterv1.h \
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
  testexternals_calo_util \
  testexternals_calo_io


testexternals_calo_util_SOURCES = testexternals.C
testexternals_calo_util_LDADD = libcalo_util.la

testexternals_calo_io_SOURCES = testexternals.C
testexternals_calo_io_LDADD = libcalo_io.la

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
