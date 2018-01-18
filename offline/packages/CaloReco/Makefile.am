AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcalo_reco.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include 
  
libcalo_reco_la_SOURCES = \
  RawClusterBuilderGraph.cc \
  RawClusterBuilderGraph_Dict.cc \
  RawClusterBuilderTemplate.cc \
  RawClusterBuilderTemplate_Dict.cc \
  RawClusterBuilderFwd.cc \
  RawClusterBuilderFwd_Dict.cc \
  RawClusterPositionCorrection.cc \
  RawClusterPositionCorrection_Dict.cc \
  RawTowerCombiner.cc \
  RawTowerCombiner_Dict.cc \
  RawTowerCalibration.cc \
  RawTowerCalibration_Dict.cc \
  BEmcCluster.cc \
  BEmcRec.cc

libcalo_reco_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcalo_reco_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lgsl \
  -lgslcblas \
  -lcalo_io \
  -lcalo_util  \
  -lphparameter
  
# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

# pkginclude_HEADERS = 

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_calo_reco


testexternals_calo_reco_SOURCES = testexternals.C
testexternals_calo_reco_LDADD = libcalo_reco.la

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
