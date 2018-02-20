AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

libPrototype2_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

lib_LTLIBRARIES = \
  libPrototype2.la

pkginclude_HEADERS = \
  RawTower_Prototype2.h \
  RawTower_Temperature.h \
  PROTOTYPE2_FEM.h

libPrototype2_la_SOURCES = \
  RawTower_Prototype2.cc \
  RawTower_Prototype2Dict.C \
  Prototype2DSTReader.cc \
  Prototype2DSTReaderDict.C \
  CaloUnpackPRDF.C \
  CaloUnpackPRDFDict.C \
  TempInfoUnpackPRDF.C \
  TempInfoUnpackPRDFDict.C \
  RunInfoUnpackPRDF.C \
  RunInfoUnpackPRDFDict.C \
  CaloCalibration.C \
  CaloCalibrationDict.C \
  GenericUnpackPRDF.C \
  GenericUnpackPRDFDict.C \
  PROTOTYPE2_FEM.C \
  PROTOTYPE2_FEMDict.C \
  RawTower_Temperature.cc \
  RawTower_TemperatureDict.C

libPrototype2_la_LIBADD = \
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
  libPrototype2.la

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
