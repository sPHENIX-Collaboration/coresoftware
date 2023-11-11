##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcalo_reco.la

#
libcalo_reco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

if USE_ONLINE

else
libcalo_reco_la_LIBADD = \
  -lmbd_io \
  -lcalo_io \
  -lcdbobjects \
  -lCLHEP \
  -lffamodules \
  -lgsl \
  -lgslcblas \
  -lglobalvertex_io \
  -lsph_onnx \
  -lphparameter \
  -lphool \
  -lSubsysReco
endif

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

if USE_ONLINE
pkginclude_HEADERS = \
  CaloWaveformFitting.h

else
pkginclude_HEADERS = \
  CaloGeomMapping.h \
  CaloWaveformFitting.h \
  CaloWaveformProcessing.h \
  CaloRecoUtility.h \
  CaloTowerBuilder.h \
  CaloTowerCalib.h \
  CaloTowerStatus.h \
  CaloTowerDefs.h \
  RawClusterBuilderGraph.h \
  RawClusterBuilderTopo.h \
  RawClusterBuilderTemplate.h \
  RawClusterPositionCorrection.h \
  RawClusterZVertexRecorrect.h \
  RawTowerCalibration.h \
  RawTowerDeadTowerInterp.h \
  DeadHotMapLoader.h \
  RawClusterDeadHotMask.h \
  TowerInfoDeadHotMask.h
endif

if USE_ONLINE
libcalo_reco_la_SOURCES = \
  CaloWaveformFitting.cc

else
libcalo_reco_la_SOURCES = \
  BEmcCluster.cc \
  BEmcProfile.cc \
  BEmcRec.cc \
  BEmcRecCEMC.cc \
  CaloGeomMapping.cc \
  CaloRecoUtility.cc \
  CaloWaveformFitting.cc \
  CaloWaveformProcessing.cc \
  CaloTowerBuilder.cc \
  CaloTowerCalib.cc \
  CaloTowerStatus.cc \
  RawClusterBuilderGraph.cc \
  RawClusterBuilderTopo.cc \
  RawClusterBuilderTemplate.cc \
  RawClusterPositionCorrection.cc \
  RawClusterZVertexRecorrect.cc \
  RawTowerCombiner.cc \
  RawTowerCalibration.cc \
  RawTowerDeadTowerInterp.cc \
  DeadHotMapLoader.cc \
  RawClusterDeadHotMask.cc \
  TowerInfoDeadHotMask.cc
endif

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_calo_reco

BUILT_SOURCES  = testexternals.cc

testexternals_calo_reco_SOURCES = testexternals.cc
testexternals_calo_reco_LDADD = libcalo_reco.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
