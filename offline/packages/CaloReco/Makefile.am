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
  -L$(OFFLINE_MAIN)/lib

libcalo_reco_la_LIBADD = \
  -lCLHEP \
  -lphool \
  -lSubsysReco \
  -lgsl \
  -lgslcblas \
  -lg4vertex_io \
  -lcalo_io \
  -lsph_onnx \
  -lcaloCalibDBFile \
  -lcdbobjects \
  -lphparameter

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  CaloWaveformProcessing.h \
  CaloTowerBuilder.h \
  CaloTowerCalib.h \
  RawClusterBuilderFwd.h \
  RawClusterBuilderGraph.h \
  RawClusterBuilderTopo.h \
  RawClusterBuilderTemplate.h \
  RawClusterPositionCorrection.h \
  RawTowerCalibration.h 

libcalo_reco_la_SOURCES = \
  BEmcCluster.cc \
  BEmcProfile.cc \
  BEmcRec.cc \
  BEmcRecCEMC.cc \
  BEmcRecFEMC.cc \
  BEmcRecEEMC.cc \
  CaloWaveformProcessing.cc \
  CaloTowerBuilder.cc \
  CaloTowerCalib.cc \
  RawClusterBuilderGraph.cc \
  RawClusterBuilderTopo.cc \
  RawClusterBuilderTemplate.cc \
  RawClusterBuilderFwd.cc \
  RawClusterPositionCorrection.cc \
  RawClusterDeadAreaMask.cc \
  RawTowerCombiner.cc \
  RawTowerCalibration.cc \
  RawTowerDeadTowerInterp.cc \
  RawTowerDeadMapLoader.cc

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
