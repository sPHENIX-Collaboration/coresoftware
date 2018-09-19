AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libtpcdaq.la

AM_CXXFLAGS = -Wall -Werror -msse2

libtpcdaq_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libtpcdaq_la_LIBADD = \
  -lfun4all \
  -lg4dst 

libtpcdaq_la_SOURCES = \
  TPCIntegratedCharge.cc \
  TPCIntegratedCharge_Dict.C \
  TPCDataStreamEmulator.cc \
  TPCDataStreamEmulator_Dict.C  \
  TPCFEETestRecov1.cc \
  TPCFEETestRecov1_Dict.C  \
  TPCDaqDefs.cc \
  TPCDaqDefs_Dict.C 

# Rule for generating table CINT dictionaries.
%_Dict.C: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  TPCDaqDefs.h 
  
################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libtpcdaq.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* testexternals.C
