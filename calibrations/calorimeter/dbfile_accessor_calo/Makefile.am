##############################################
# please add new classes in alphabetical order

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include


AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcaloCalibDBFile.la

libcaloCalibDBFile_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libcaloCalibDBFile_la_LIBADD = \
  -lphool 


pkginclude_HEADERS = \
  CaloCalibSimpleCorrFile.h \
  HcalCaloCalibSimpleCorrFilev1.h \
  CEmcCaloCalibSimpleCorrFilev1.h


libcaloCalibDBFile_la_SOURCES = \
  HcalCaloCalibSimpleCorrFilev1.cc \
  CEmcCaloCalibSimpleCorrFilev1.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libcaloCalibDBFile.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
