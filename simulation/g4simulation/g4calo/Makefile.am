AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libg4calo.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include


libg4calo_la_SOURCES = \
  HcalRawTowerBuilder.cc \
  HcalRawTowerBuilder_Dict.cc \
  Prototype2RawTowerBuilder.cc \
  Prototype2RawTowerBuilder_Dict.cc \
  RawTowerBuilder.cc \
  RawTowerBuilder_Dict.cc \
  RawTowerBuilderByHitIndex.cc \
  RawTowerBuilderByHitIndex_Dict.cc \
  RawTowerDigitizer.cc \
  RawTowerDigitizer_Dict.cc

libg4calo_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libg4calo_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lg4detectors_io \
  -lg4testbench \
  -lphg4hit \
  -lphparameter \
  -lgsl \
  -lgslcblas \
  -lcalo_io \
  -lcalo_util

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  HcalRawTowerBuilder.h \
  RawTowerBuilder.h \
  RawTowerDigitizer.h

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4calo


testexternals_g4calo_SOURCES = testexternals.C
testexternals_g4calo_LDADD = libg4calo.la

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
