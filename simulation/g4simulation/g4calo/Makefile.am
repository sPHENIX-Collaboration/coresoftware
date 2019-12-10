##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libg4calo.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I${G4_MAIN}/include

pkginclude_HEADERS = \
  HcalRawTowerBuilder.h \
  RawTowerBuilder.h \
  RawTowerBuilderByHitIndex.h \
  RawTowerDigitizer.h

if ! MAKEROOT6
  ROOT5_DICTS = \
    HcalRawTowerBuilder_Dict.cc \
    RawTowerBuilder_Dict.cc \
    RawTowerBuilderByHitIndex_Dict.cc \
    RawTowerDigitizer_Dict.cc
endif

libg4calo_la_SOURCES = \
  $(ROOT5_DICTS) \
  HcalRawTowerBuilder.cc \
  RawTowerBuilder.cc \
  RawTowerBuilderByHitIndex.cc \
  RawTowerDigitizer.cc

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

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_g4calo

BUILT_SOURCES = testexternals.cc

testexternals_g4calo_SOURCES = testexternals.cc
testexternals_g4calo_LDADD = libg4calo.la

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
