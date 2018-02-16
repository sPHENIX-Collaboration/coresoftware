AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

libPrototype3_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

lib_LTLIBRARIES = \
  libPrototype3.la

pkginclude_HEADERS = \
  RawTower_Prototype3.h \
  RawTower_Temperature.h \
  PROTOTYPE3_FEM.h

libPrototype3_la_SOURCES = \
  RawTower_Prototype3.cc \
  RawTower_Prototype3Dict.C \
  Prototype3DSTReader.cc \
  Prototype3DSTReaderDict.C \
  CaloUnpackPRDF.C \
  CaloUnpackPRDFDict.C \
  TempInfoUnpackPRDF.C \
  TempInfoUnpackPRDFDict.C \
  RunInfoUnpackPRDF.C \
  RunInfoUnpackPRDFDict.C \
  EventInfoSummary.C \
  EventInfoSummaryDict.C \
  CaloCalibration.C \
  CaloCalibrationDict.C \
  GenericUnpackPRDF.C \
  GenericUnpackPRDFDict.C \
  PROTOTYPE3_FEM.C \
  PROTOTYPE3_FEMDict.C \
  RawTower_Temperature.cc \
  RawTower_TemperatureDict.C

libPrototype3_la_LIBADD = \
  -lSubsysReco \
  -lphool \
  -lfun4all \
  -lcalo_util \
  -lphparameter

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPrototype3.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order
# Rule for generating CINT dictionaries from class headers.
%Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *Dict*

testexternals_SOURCES = testexternals.C
