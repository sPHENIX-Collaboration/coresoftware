AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS=-std=c++0x 

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcemc.la

INCLUDES = \
  -I$(includedir) \
  -isystem $(OFFLINE_MAIN)/include \
  -isystem $(ROOTSYS)/include

RINCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

libcemc_la_SOURCES = \
  Cal_Dict.C \
  RawClusterBuilder.cc \
  RawClusterBuilderv1.cc \
  RawClusterv1.cc \
  RawClusterContainer.cc \
  RawTower.cc \
  RawTowerv1.cc \
  RawTowerv2.cc \
  RawTowerBuilder.cc \
  RawTowerBuilderCone.cc \
  RawTowerCombiner.cc \
  RawTowerContainer.cc \
  RawTowerGeom.cc \
  RawTowerGeomv1.cc \
  RawTowerGeomv2.cc \
  TowerDigitization.cc \
  BEmcRec.cc \
  BEmcCluster.cc


libcemc_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcemc_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lg4detectors \
  -lgsl \
  -lgslcblas

Cal_Dict.C: \
  RawClusterBuilder.h \
  RawClusterBuilderv1.h \
  RawCluster.h \
  RawClusterv1.h \
  RawClusterContainer.h \
  RawTower.h \
  RawTowerv1.h \
  RawTowerv2.h \
  RawTowerBuilder.h \
  RawTowerBuilderCone.h \
  RawTowerCombiner.h \
  RawTowerContainer.h \
  RawTowerGeom.h \
  RawTowerGeomv1.h \
  RawTowerGeomv2.h \
  TowerDigitization.h \
  BEmcRec.h \
  BEmcCluster.h \
  CalLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(RINCLUDES) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  RawCluster.h \
  RawClusterContainer.h \
  RawTower.h \
  RawTowerv1.h \
  RawTowerv2.h \
  RawTowerContainer.h \
  RawTowerGeom.h


################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_cemc 

testexternals_cemc_SOURCES = testexternals.C
testexternals_cemc_LDADD = libcemc.la


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
